"""v46 = v45 with Bochner-Toeplitz extended to degree 2K.

Adds auxiliary Fourier-coefficient variables a[K+1..2K], b[K+1..2K]
with their own cell bounds, and replaces v45's (K+1)-Bochner with
a (2K+1)-Bochner that ties all of a[0..2K], b[0..2K] together.

The Fejer-Riesz constraint still uses only a[1..K], b[1..K] (via M).
"""

from __future__ import annotations

from dataclasses import dataclass

import cvxpy as cp
import numpy as np


@dataclass
class V46Result:
    status: str
    Omega: float
    primal: float


class AutocorrLowerBoundV46:
    def __init__(self, N: int = 8, K: int = 8) -> None:
        self.N = N
        self.K = K
        self.K2 = 2 * K  # extended order for Bochner-Toeplitz
        self.dim_p = 2 * N
        self.L_f = 1.0 / (4 * N)

        self.Omega = cp.Variable(nonneg=True, name="Omega")
        self.p = cp.Variable(self.dim_p, nonneg=True, name="p")
        # Extended a, b up to 2K
        self.a = cp.Variable(self.K2 + 1, name="a")
        self.b = cp.Variable(self.K2 + 1, name="b")
        # Moment matrix M still on (a_1..K, b_1..K)
        self.M = cp.Variable((2 * K + 1, 2 * K + 1), symmetric=True, name="M")
        self.v = cp.Variable(K, nonneg=True, name="v")
        # Fejer-Riesz Q for sigma_K
        self.Q = cp.Variable((K + 1, K + 1), hermitian=True, name="Q")

        constraints: list = []

        # Mass normalization
        constraints.append(cp.sum(self.p) == 1)
        constraints.append(self.a[0] == 1)
        constraints.append(self.b[0] == 0)

        # Moment lift M (only k=0..K)
        constraints.append(self.M >> 0)
        constraints.append(self.M[0, 0] == 1)
        for k in range(1, K + 1):
            constraints.append(self.M[0, k] == self.a[k])
            constraints.append(self.M[k, 0] == self.a[k])
            constraints.append(self.M[0, K + k] == self.b[k])
            constraints.append(self.M[K + k, 0] == self.b[k])

        # Cell bounds on a_k, b_k for k=1..2K
        cos_min = np.zeros((self.K2 + 1, self.dim_p))
        cos_max = np.zeros((self.K2 + 1, self.dim_p))
        sin_min = np.zeros((self.K2 + 1, self.dim_p))
        sin_max = np.zeros((self.K2 + 1, self.dim_p))
        for k in range(1, self.K2 + 1):
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
        for k in range(1, self.K2 + 1):
            constraints.append(self.a[k] >= cos_min[k] @ self.p)
            constraints.append(self.a[k] <= cos_max[k] @ self.p)
            constraints.append(self.b[k] >= sin_min[k] @ self.p)
            constraints.append(self.b[k] <= sin_max[k] @ self.p)

        # v14 squared autocorrelation Fejer (k=1..K only, redundant)
        for k in range(1, K + 1):
            sum_diag = self.M[k, k] + self.M[K + k, K + k]
            constraints.append(self.v[k - 1] >= cp.square(sum_diag))
        weights = np.array([1.0 - k / (K + 1) for k in range(1, K + 1)])
        constraints.append(self.Omega >= 1.0 + 2.0 * (weights @ self.v))

        # v44 Fejer-Riesz exact (degree K)
        constraints.append(self.Q >> 0)
        constraints.append(cp.real(cp.trace(self.Q)) == self.Omega - 1)
        for k in range(1, K + 1):
            anti_diag_sum = sum(self.Q[i, i + k] for i in range(K + 1 - k))
            constraints.append(
                cp.real(anti_diag_sum)
                == -weights[k - 1] * (self.M[k, k] - self.M[K + k, K + k])
            )
            constraints.append(
                cp.imag(anti_diag_sum) == 2.0 * weights[k - 1] * self.M[k, K + k]
            )

        # NEW v46: extended Bochner-Toeplitz of size (2K+1)x(2K+1)
        T_R_rows = []
        T_I_rows = []
        for i in range(self.K2 + 1):
            row_R = []
            row_I = []
            for j in range(self.K2 + 1):
                d = abs(i - j)
                if d == 0:
                    row_R.append(self.a[0])
                    row_I.append(cp.Constant(0.0))
                else:
                    row_R.append(self.a[d])
                    if i > j:
                        row_I.append(-self.b[d])
                    else:
                        row_I.append(self.b[d])
            T_R_rows.append(cp.hstack(row_R))
            T_I_rows.append(cp.hstack(row_I))
        T_R = cp.vstack(T_R_rows)
        T_I = cp.vstack(T_I_rows)
        big = cp.bmat([[T_R, -T_I], [T_I, T_R]])
        constraints.append(big >> 0)

        self.constraints = constraints
        self.problem = cp.Problem(cp.Minimize(self.Omega), constraints)

    def solve(self, solver: str = "MOSEK", verbose: bool = False, **kwargs) -> V46Result:
        val = self.problem.solve(solver=solver, verbose=verbose, **kwargs)
        return V46Result(
            status=self.problem.status,
            Omega=float(self.Omega.value) if self.Omega.value is not None else float("nan"),
            primal=float(val) if val is not None else float("nan"),
        )


if __name__ == "__main__":
    import time

    print("v46 (extended Bochner-Toeplitz to 2K) sanity sweep:")
    for N, K in [(8, 8), (16, 16), (16, 32), (24, 24), (32, 32), (50, 32)]:
        t = time.time()
        try:
            prob = AutocorrLowerBoundV46(N=N, K=K)
            out = prob.solve(solver="MOSEK", verbose=False)
            print(f"N={N:3d}, K={K:3d}: Omega={out.Omega:.6f}  ({time.time()-t:.1f}s)")
        except Exception as e:
            print(f"N={N}, K={K}: failed: {e}")
