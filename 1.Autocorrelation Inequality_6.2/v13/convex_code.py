"""v13: v12 + Bochner / Toeplitz PSD constraint on f's Fourier sequence.

Bochner's theorem says that the sequence (f_hat(k))_k is positive
definite (i.e., the Hermitian Toeplitz matrix [f_hat(j-k)]_{j,k>=0} is
PSD) if and only if f >= 0 as a (Borel) measure.  Since admissible f
is nonneg, this constraint must hold.

In real form, with f_hat(k) = a_k - i b_k (k >= 1) and f_hat(0) = 1,
f_hat(-k) = a_k + i b_k, the Hermitian Toeplitz is

   T = T_R + i T_I,    T_R[i,j] = a_{|i-j|},  T_I[i,j] = sign(j-i)*b_{|i-j|}.

The Hermitian PSD constraint is equivalent to the real PSD constraint
on the (2(K+1)) x (2(K+1)) block matrix
   [[T_R, -T_I],
    [T_I,  T_R]]   >> 0.

We add this to the v12 SDP.

The expectation: this LMI on (a, b) constrains the lifted moments
M[ic,ic], M[K+k,K+k], M[ic,id] indirectly (through PSD-on-(p,p^T)
joint with the lifted M).
"""

from __future__ import annotations

from dataclasses import dataclass

import cvxpy as cp
import numpy as np


@dataclass
class V13Result:
    status: str
    Omega: float
    primal: float


class AutocorrLowerBoundV13:
    def __init__(self, N: int = 8, K: int = 6) -> None:
        if N < 1 or K < 1:
            raise ValueError("N, K must be positive")
        self.N = N
        self.K = K
        self.dim_p = 2 * N
        self.dim_W = 2 * N
        self.L_f = 1.0 / (4 * N)
        self.L_F = 1.0 / (2 * N)

        self.Omega = cp.Variable(nonneg=True, name="Omega")
        self.p = cp.Variable(self.dim_p, nonneg=True, name="p")
        self.W = cp.Variable(self.dim_W, nonneg=True, name="W")
        self.a = cp.Variable(K + 1, name="a")
        self.b = cp.Variable(K + 1, name="b")
        self.M = cp.Variable((2 * K + 1, 2 * K + 1), symmetric=True, name="M")

        constraints: list = []

        # mass and basics
        constraints.append(cp.sum(self.p) == 1)
        constraints.append(self.L_F * cp.sum(self.W) == 1)
        constraints.append(self.W <= self.Omega)
        constraints.append(self.a[0] == 1)
        constraints.append(self.b[0] == 0)

        # PSD lift on M
        constraints.append(self.M >> 0)
        constraints.append(self.M[0, 0] == 1)
        for k in range(1, K + 1):
            constraints.append(self.M[0, k] == self.a[k])
            constraints.append(self.M[k, 0] == self.a[k])
            constraints.append(self.M[0, K + k] == self.b[k])
            constraints.append(self.M[K + k, 0] == self.b[k])

        # f-side bounds on a_k, b_k from p
        cos_min_f = np.zeros((K + 1, self.dim_p))
        cos_max_f = np.zeros((K + 1, self.dim_p))
        sin_min_f = np.zeros((K + 1, self.dim_p))
        sin_max_f = np.zeros((K + 1, self.dim_p))
        for k in range(1, K + 1):
            for j in range(self.dim_p):
                a_l = -0.25 + j * self.L_f
                a_r = a_l + self.L_f
                xs = np.linspace(a_l, a_r, 401)
                cv = np.cos(2 * np.pi * k * xs)
                sv = np.sin(2 * np.pi * k * xs)
                cos_min_f[k, j] = cv.min()
                cos_max_f[k, j] = cv.max()
                sin_min_f[k, j] = sv.min()
                sin_max_f[k, j] = sv.max()
        for k in range(1, K + 1):
            constraints.append(self.a[k] >= cos_min_f[k] @ self.p)
            constraints.append(self.a[k] <= cos_max_f[k] @ self.p)
            constraints.append(self.b[k] >= sin_min_f[k] @ self.p)
            constraints.append(self.b[k] <= sin_max_f[k] @ self.p)

        # F-side bounds on F_R(k), F_I(k) from W
        cos_min_F = np.zeros((K + 1, self.dim_W))
        cos_max_F = np.zeros((K + 1, self.dim_W))
        sin_min_F = np.zeros((K + 1, self.dim_W))
        sin_max_F = np.zeros((K + 1, self.dim_W))
        for k in range(1, K + 1):
            for m in range(self.dim_W):
                t_l = -0.5 + m * self.L_F
                t_r = t_l + self.L_F
                ts = np.linspace(t_l, t_r, 401)
                cv = np.cos(2 * np.pi * k * ts)
                sv = np.sin(2 * np.pi * k * ts)
                cos_min_F[k, m] = cv.min()
                cos_max_F[k, m] = cv.max()
                sin_min_F[k, m] = sv.min()
                sin_max_F[k, m] = sv.max()
        for k in range(1, K + 1):
            F_R_k = self.M[k, k] - self.M[K + k, K + k]
            F_I_k = -2.0 * self.M[k, K + k]
            constraints.append(F_R_k >= self.L_F * (cos_min_F[k] @ self.W))
            constraints.append(F_R_k <= self.L_F * (cos_max_F[k] @ self.W))
            constraints.append(-F_I_k >= self.L_F * (sin_min_F[k] @ self.W))
            constraints.append(-F_I_k <= self.L_F * (sin_max_F[k] @ self.W))

        # NEW v13: Bochner / Toeplitz PSD constraint on (a_k, b_k)
        # Build T_R and T_I as (K+1)x(K+1) matrices.
        T_R_rows = []
        T_I_rows = []
        for i in range(K + 1):
            row_R = []
            row_I = []
            for j in range(K + 1):
                d = abs(i - j)
                if d == 0:
                    row_R.append(self.a[0])
                    row_I.append(cp.Constant(0.0))
                else:
                    row_R.append(self.a[d])
                    sign = -1.0 if i > j else 1.0
                    row_I.append(sign * self.b[d])
            T_R_rows.append(cp.hstack(row_R))
            T_I_rows.append(cp.hstack(row_I))
        T_R = cp.vstack(T_R_rows)
        T_I = cp.vstack(T_I_rows)
        # Hermitian PSD as real LMI: [[T_R, -T_I], [T_I, T_R]] >> 0
        big = cp.bmat([[T_R, -T_I], [T_I, T_R]])
        constraints.append(big >> 0)

        self.constraints = constraints
        self.problem = cp.Problem(cp.Minimize(self.Omega), constraints)

    def solve(self, solver: str = "MOSEK", verbose: bool = False, **kwargs) -> V13Result:
        val = self.problem.solve(solver=solver, verbose=verbose, **kwargs)
        return V13Result(
            status=self.problem.status,
            Omega=float(self.Omega.value) if self.Omega.value is not None else float("nan"),
            primal=float(val) if val is not None else float("nan"),
        )


if __name__ == "__main__":
    for N in (8, 12):
        for K in (6, 12):
            try:
                prob = AutocorrLowerBoundV13(N=N, K=K)
                out = prob.solve(solver="MOSEK", verbose=False)
                print(f"N={N:2d}, K={K:2d}: Omega={out.Omega:.6f}")
            except Exception as e:
                print(f"N={N}, K={K}: failed: {e}")
