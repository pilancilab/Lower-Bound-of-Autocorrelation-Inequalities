"""v12: F-side LP with Bochner Toeplitz on F's Fourier sequence.

Per the v7 verdict feedback item 2:
"Don't restrict to even f: drop (E)/(F)/(G), and find a different way
to break the diagonal-Q degeneracy.  The most promising is the F-side
autocorrelation Toeplitz lift.  Define cell averages W_m of F = f*f
directly on [-1/2, 1/2], use the rank-one identification
F_hat(k) = f_hat(k)^2, and add the Toeplitz PSD constraint on the
discrete autocorrelation F_hat. This is a UNIVERSAL constraint
(every admissible f yields a Toeplitz-PSD autocorrelation), not
restricted to even f."

The idea
--------

Variables (general, not even):
* Omega.
* W_m: cell averages of F = f*f on cells of [-1/2, 1/2].
* a_k, b_k: Fourier coefficients of f (k = 0..K), with a_0 = 1, b_0 = 0.
* M: PSD lift on (1, a_1..a_K, b_1..b_K) of size 2K+1.
* F_R, F_I: real/imaginary parts of F_hat(k), k = -K..K (real F gives
  F_I(-k) = -F_I(k), F_R(-k) = F_R(k); we store only k = 0..K).

Constraints:
1. W_m >= 0, sum_m W_m * L_F = 1, W_m <= Omega.
2. Linear bounds on F_R(k), F_I(k) from W_m via interval bounds:
   F_R(k) = int F * cos(2*pi*k*t) dt is sandwiched between
   sum_m W_m * (min cos)_{m,k} * L_F  and  sum_m W_m * (max cos)_{m,k} * L_F.
   Similarly F_I(k) = -int F * sin(2*pi*k*t) dt.
3. Rank-one identification F_hat(k) = f_hat(k)^2:
   F_R(k) = (a_k^2 - b_k^2)  -- lifted as M[k,k] - M[K+k, K+k]
   F_I(k) = -2 a_k b_k       -- lifted as -2 * M[k, K+k]
4. PSD lift on M: M >> 0, M[0,0]=1, M[0, k]=a_k, M[0, K+k]=b_k.
5. Bochner / Toeplitz PSD on F_hat: [F_hat(j-k)]_{j,k=0..K} >> 0
   (Hermitian PSD, expressed as a real block matrix LMI).

The Bochner constraint encodes that there is a NONNEG function F
(not just any function) consistent with the F_hat values.  It is
implied by W_m >= 0, but adds structural information when combined
with the rank-one link to f-side variables.

Status: this is an experimental implementation.  Numerical performance
is reported in rigorousproof.md.
"""

from __future__ import annotations

from dataclasses import dataclass

import cvxpy as cp
import numpy as np


@dataclass
class V12Result:
    status: str
    Omega: float
    primal: float


class AutocorrLowerBoundV12:
    def __init__(self, N: int = 8, K: int = 6) -> None:
        if N < 1 or K < 1:
            raise ValueError("N, K must be positive")
        self.N = N
        self.K = K
        self.dim_p = 2 * N
        self.dim_W = 2 * N
        self.L_f = 1.0 / (4 * N)  # f-cell width
        self.L_F = 1.0 / (2 * N)  # F-cell width

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
        # F_R(k) = int F cos(2 pi k t) dt
        # F_I(k) = -int F sin(2 pi k t) dt
        # rank-one identification: F_R(k) = a_k^2 - b_k^2 = M[k,k] - M[K+k, K+k]
        #                          F_I(k) = -2 a_k b_k = -2 * M[k, K+k]
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

        self.constraints = constraints
        self.problem = cp.Problem(cp.Minimize(self.Omega), constraints)

    def solve(self, solver: str = "MOSEK", verbose: bool = False, **kwargs) -> V12Result:
        val = self.problem.solve(solver=solver, verbose=verbose, **kwargs)
        return V12Result(
            status=self.problem.status,
            Omega=float(self.Omega.value) if self.Omega.value is not None else float("nan"),
            primal=float(val) if val is not None else float("nan"),
        )


if __name__ == "__main__":
    for N in (8, 12):
        for K in (6, 12):
            try:
                prob = AutocorrLowerBoundV12(N=N, K=K)
                out = prob.solve(solver="MOSEK", verbose=False)
                print(f"N={N:2d}, K={K:2d}: Omega={out.Omega:.6f}")
            except Exception as e:
                print(f"N={N}, K={K}: failed: {e}")
