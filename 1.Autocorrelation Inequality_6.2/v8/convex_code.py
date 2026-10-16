"""v8 convex relaxation: F-side LP with Fourier coupling.

Builds on v7's rigorous (W) constraints and adds the F-side variables
W_m (cell averages of F = f*f on [-1/2, 1/2]) plus a Fourier coupling
linking f's Fourier coefficients to F's via the squared identity.

The relaxation bounds C^{even}_{6.2} (the even-restricted constant),
addressing the verifier's recommendation to switch to the F-side.

Variables
---------

* Omega                       (the value to minimize, lower bound on
                               C^{even}_{6.2}).
* p_j   (j = 1..2N)            mass of f on the j-th interval of [-1/4, 1/4],
                               with p_j = p_{2N+1-j} (even-f restriction).
* W_m   (m = 1..2N)            avg of F = f*f on the m-th interval of [-1/2, 1/2],
                               with W_m = W_{2N+1-m} (F is even when f is even).
* a_k   (k = 1..K)             cosine Fourier coefficients of f
                               (sine coefficients vanish for even f).
* M_kk  (k = 1..K)             SDP-lifted variables for a_k^2.

Constraints
-----------

(P0) Mass: sum p = 1, sum_m W_m * L_F = 1.
(P1) Even symmetry: p_j = p_{2N+1-j}, W_m = W_{2N+1-m}.
(P2) f-side bounds on a_k from p:
       sum_j p_j * (min cos)_{j,k}  <=  a_k  <=  sum_j p_j * (max cos)_{j,k}
(P3) PSD lift on (1, a_1, ..., a_K):
       M_kk >= a_k^2  (rotated SOC)  AND  the full M >> 0 with M[0,0]=1.
(P4) Parseval lower bound:  sum_k M_kk  >=  1/2
       (this is a LINEAR constraint on M_kk; valid because for
        every admissible f, sum_k a_k^2 = (int f^2 - 1)/2 >= 1/2
        by Cauchy-Schwarz on the [-1/4, 1/4] support.)
(P5) F-side Fourier linkage:
       A_k(F) = M_kk  (for even f: A_k = a_k^2 - b_k^2 = a_k^2;
                        in lift, M_kk lifts a_k^2)
       A_k(F) is bounded by W via interval bounds:
         L_F * sum_m W_m * (min cos)_{m, k}  <=  M_kk  <=  L_F * sum_m W_m * (max cos)_{m, k}.
(P6) Omega constraint:
       W_m  <=  Omega  for all m.

The Parseval constraint (P4) is what breaks the diagonal-Q-style
degeneracy on the f-side: for any "mixture" measure on the simplex,
the average a_k^2 over the mixture obeys Jensen's inequality, so
sum_k <a_k^2> = 1/2 forces the mixture's mass to be spread out in
Fourier space, which incompatible with the F-side bound (a_k^2 = M_kk
must be small for uniform W).

Validity is established in rigorousproof.md.
"""

from __future__ import annotations

from dataclasses import dataclass

import cvxpy as cp
import numpy as np


@dataclass
class V8Result:
    status: str
    Omega: float
    primal: float


class AutocorrLowerBoundV8:
    def __init__(self, N: int = 8, K: int = 6) -> None:
        if N < 1 or K < 1:
            raise ValueError("N and K must be positive integers")
        self.N = N
        self.K = K
        self.dim_p = 2 * N
        self.dim_W = 2 * N
        self.L_f = 1.0 / (4 * N)  # f-interval width on [-1/4, 1/4]
        self.L_F = 1.0 / (2 * N)  # F-interval width on [-1/2, 1/2]

        self.Omega = cp.Variable(nonneg=True, name="Omega")
        self.p = cp.Variable(self.dim_p, nonneg=True, name="p")
        self.W = cp.Variable(self.dim_W, nonneg=True, name="W")
        self.a = cp.Variable(K + 1, name="a")  # a_0 = 1, a_1..a_K
        # M is the (K+1)x(K+1) PSD lift on (1, a_1, ..., a_K)
        self.M = cp.Variable((K + 1, K + 1), symmetric=True, name="M")

        constraints: list = []

        # ---- (P0) mass + (P6) Omega bound on W ----
        constraints.append(cp.sum(self.p) == 1)
        constraints.append(self.L_F * cp.sum(self.W) == 1)
        constraints.append(self.W <= self.Omega)
        constraints.append(self.a[0] == 1)

        # ---- (P1) even symmetry ----
        for i in range(self.dim_p // 2):
            constraints.append(self.p[i] == self.p[self.dim_p - 1 - i])
        for m in range(self.dim_W // 2):
            constraints.append(self.W[m] == self.W[self.dim_W - 1 - m])

        # ---- (P3) M PSD lift ----
        constraints.append(self.M >> 0)
        constraints.append(self.M[0, 0] == 1)
        for k in range(1, K + 1):
            constraints.append(self.M[0, k] == self.a[k])
            constraints.append(self.M[k, 0] == self.a[k])

        # ---- (P2) f-side bounds on a_k from p ----
        cos_min_f = np.zeros((K + 1, self.dim_p))
        cos_max_f = np.zeros((K + 1, self.dim_p))
        for k in range(1, K + 1):
            for j in range(self.dim_p):
                a_left = -0.25 + j * self.L_f
                a_right = a_left + self.L_f
                xs = np.linspace(a_left, a_right, 401)
                cv = np.cos(2.0 * np.pi * k * xs)
                cos_min_f[k, j] = cv.min()
                cos_max_f[k, j] = cv.max()
        for k in range(1, K + 1):
            constraints.append(self.a[k] >= cos_min_f[k] @ self.p)
            constraints.append(self.a[k] <= cos_max_f[k] @ self.p)

        # ---- (P4) Parseval constraint ----
        # sum_{k=1}^K M[k,k] >= 1/2
        # (This is sum a_k^2 >= 1/2 lifted; valid for any admissible f.)
        constraints.append(
            cp.sum(cp.hstack([self.M[k, k] for k in range(1, K + 1)])) >= 0.5
        )

        # ---- (P5) F-side bound on A_k from W, with A_k = M[k,k] (even f) ----
        cos_min_F = np.zeros((K + 1, self.dim_W))
        cos_max_F = np.zeros((K + 1, self.dim_W))
        for k in range(1, K + 1):
            for m in range(self.dim_W):
                t_left = -0.5 + m * self.L_F
                t_right = t_left + self.L_F
                ts = np.linspace(t_left, t_right, 401)
                cv = np.cos(2.0 * np.pi * k * ts)
                cos_min_F[k, m] = cv.min()
                cos_max_F[k, m] = cv.max()
        for k in range(1, K + 1):
            constraints.append(self.M[k, k] >= self.L_F * (cos_min_F[k] @ self.W))
            constraints.append(self.M[k, k] <= self.L_F * (cos_max_F[k] @ self.W))

        self.constraints = constraints
        self.problem = cp.Problem(cp.Minimize(self.Omega), constraints)

    def solve(self, solver: str = "MOSEK", verbose: bool = False, **kwargs) -> V8Result:
        val = self.problem.solve(solver=solver, verbose=verbose, **kwargs)
        return V8Result(
            status=self.problem.status,
            Omega=float(self.Omega.value) if self.Omega.value is not None else float("nan"),
            primal=float(val) if val is not None else float("nan"),
        )


if __name__ == "__main__":
    for N in (4, 6, 8):
        for K in (4, 6, 8):
            prob = AutocorrLowerBoundV8(N=N, K=K)
            try:
                out = prob.solve(solver="MOSEK", verbose=False)
            except Exception as exc:
                print(f"N={N}, K={K}: MOSEK failed: {exc}")
                continue
            print(f"N={N:2d}, K={K:2d}: status={out.status}, Omega={out.Omega:.6f}")
