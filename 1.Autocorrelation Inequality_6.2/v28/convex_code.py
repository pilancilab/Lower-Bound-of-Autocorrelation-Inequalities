"""v28: v14 + Bochner Toeplitz on F's Fourier sequence (NOT f's).

Bochner: F >= 0 implies the Hermitian Toeplitz matrix
[F_hat(j-k)]_{j,k=0..K_F} >> 0.  In our lift, F_hat(k) = (a_k - i b_k)^2,
so its real and imaginary parts are
   Re F_hat(k) = a_k^2 - b_k^2 = M[k,k] - M[K+k, K+k]
   Im F_hat(k) = -2 a_k b_k = -2 M[k, K+k]
both LINEAR in the lifted M.

The Hermitian PSD constraint is expressed as a real block LMI.
"""

from __future__ import annotations

from dataclasses import dataclass

import cvxpy as cp
import numpy as np


@dataclass
class V28Result:
    status: str
    Omega: float
    primal: float


class AutocorrLowerBoundV28:
    def __init__(self, N: int = 8, K: int = 8) -> None:
        if N < 1 or K < 1:
            raise ValueError("N, K must be positive")
        self.N = N
        self.K = K
        self.dim_p = 2 * N
        self.L_f = 1.0 / (4 * N)

        self.Omega = cp.Variable(nonneg=True, name="Omega")
        self.p = cp.Variable(self.dim_p, nonneg=True, name="p")
        self.a = cp.Variable(K + 1, name="a")
        self.b = cp.Variable(K + 1, name="b")
        self.M = cp.Variable((2 * K + 1, 2 * K + 1), symmetric=True, name="M")
        self.v = cp.Variable(K, nonneg=True, name="v")

        constraints: list = []

        constraints.append(cp.sum(self.p) == 1)
        constraints.append(self.a[0] == 1)
        constraints.append(self.b[0] == 0)

        constraints.append(self.M >> 0)
        constraints.append(self.M[0, 0] == 1)
        for k in range(1, K + 1):
            constraints.append(self.M[0, k] == self.a[k])
            constraints.append(self.M[k, 0] == self.a[k])
            constraints.append(self.M[0, K + k] == self.b[k])
            constraints.append(self.M[K + k, 0] == self.b[k])

        cos_min = np.zeros((K + 1, self.dim_p))
        cos_max = np.zeros((K + 1, self.dim_p))
        sin_min = np.zeros((K + 1, self.dim_p))
        sin_max = np.zeros((K + 1, self.dim_p))
        for k in range(1, K + 1):
            for j in range(self.dim_p):
                a_l = -0.25 + j * self.L_f
                a_r = a_l + self.L_f
                xs = np.linspace(a_l, a_r, 401)
                cv = np.cos(2 * np.pi * k * xs)
                sv = np.sin(2 * np.pi * k * xs)
                cos_min[k, j] = cv.min()
                cos_max[k, j] = cv.max()
                sin_min[k, j] = sv.min()
                sin_max[k, j] = sv.max()
        for k in range(1, K + 1):
            constraints.append(self.a[k] >= cos_min[k] @ self.p)
            constraints.append(self.a[k] <= cos_max[k] @ self.p)
            constraints.append(self.b[k] >= sin_min[k] @ self.p)
            constraints.append(self.b[k] <= sin_max[k] @ self.p)

        for k in range(1, K + 1):
            sum_diag = self.M[k, k] + self.M[K + k, K + k]
            constraints.append(self.v[k - 1] >= cp.square(sum_diag))

        weights = np.array([1.0 - k / (K + 1) for k in range(1, K + 1)])
        constraints.append(self.Omega >= 1.0 + 2.0 * (weights @ self.v))

        # NEW v28: Bochner Toeplitz on F's Fourier sequence
        # F_hat(0) = 1, F_hat(k) = (a_k^2 - b_k^2) - 2i a_k b_k for k > 0.
        # Lifted:  Re F_hat(k) = M[k,k] - M[K+k, K+k]
        #          Im F_hat(k) = -2 M[k, K+k]
        # F_hat(-k) = conj(F_hat(k)).
        # Toeplitz matrix T[i,j] = F_hat(j - i) (Hermitian convention).
        # T_R[i,j] = Re F_hat(j-i),  T_I[i,j] = Im F_hat(j-i)
        # T = T_R + i T_I, Hermitian, PSD.
        K_T = K  # Toeplitz size = K_T + 1
        T_R_rows = []
        T_I_rows = []
        for i in range(K_T + 1):
            row_R = []
            row_I = []
            for j in range(K_T + 1):
                d = j - i
                if d == 0:
                    row_R.append(cp.Constant(1.0))  # F_hat(0) = 1
                    row_I.append(cp.Constant(0.0))
                else:
                    abs_d = abs(d)
                    if abs_d > K:
                        row_R.append(cp.Constant(0.0))
                        row_I.append(cp.Constant(0.0))
                    else:
                        re_F = self.M[abs_d, abs_d] - self.M[K + abs_d, K + abs_d]
                        im_F = -2.0 * self.M[abs_d, K + abs_d]
                        row_R.append(re_F)
                        # F_hat(-k) = conj F_hat(k), so Im flips sign
                        if d < 0:
                            row_I.append(-im_F)
                        else:
                            row_I.append(im_F)
            T_R_rows.append(cp.hstack(row_R))
            T_I_rows.append(cp.hstack(row_I))
        T_R = cp.vstack(T_R_rows)
        T_I = cp.vstack(T_I_rows)
        # Real LMI: [[T_R, -T_I], [T_I, T_R]] >> 0
        big = cp.bmat([[T_R, -T_I], [T_I, T_R]])
        constraints.append(big >> 0)

        self.constraints = constraints
        self.problem = cp.Problem(cp.Minimize(self.Omega), constraints)

    def solve(self, solver: str = "MOSEK", verbose: bool = False, **kwargs) -> V28Result:
        val = self.problem.solve(solver=solver, verbose=verbose, **kwargs)
        return V28Result(
            status=self.problem.status,
            Omega=float(self.Omega.value) if self.Omega.value is not None else float("nan"),
            primal=float(val) if val is not None else float("nan"),
        )


if __name__ == "__main__":
    for N in (8, 16):
        for K in (8, 16):
            try:
                prob = AutocorrLowerBoundV28(N=N, K=K)
                out = prob.solve(solver="MOSEK", verbose=False)
                print(f"N={N:2d}, K={K:2d}: Omega={out.Omega:.6f}")
            except Exception as e:
                print(f"N={N}, K={K}: failed: {e}")
