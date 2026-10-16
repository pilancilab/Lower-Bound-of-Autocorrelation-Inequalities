"""v40: v30 + v13 Bochner Toeplitz on f's Fourier."""

from __future__ import annotations

from dataclasses import dataclass

import cvxpy as cp
import numpy as np


@dataclass
class V40Result:
    status: str
    Omega: float
    primal: float


class AutocorrLowerBoundV40:
    def __init__(self, N: int = 16, K: int = 16, n_tests: int = 51) -> None:
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

        # v30 pointwise Fejer
        ts = np.linspace(-0.5, 0.5, n_tests)
        for t_star in ts:
            cos_vals = np.cos(2 * np.pi * np.arange(1, K + 1) * t_star)
            sin_vals = np.sin(2 * np.pi * np.arange(1, K + 1) * t_star)
            terms = []
            for k in range(1, K + 1):
                re_F = self.M[k, k] - self.M[K + k, K + k]
                neg_im_F = 2.0 * self.M[k, K + k]
                weight_k = 1.0 - k / (K + 1)
                terms.append(weight_k * (re_F * cos_vals[k - 1] + neg_im_F * sin_vals[k - 1]))
            sigma_K_t = 1.0 + 2.0 * cp.sum(cp.hstack(terms))
            constraints.append(self.Omega >= sigma_K_t)

        # v13 Bochner Toeplitz on (a, b) for f's Fourier
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
        big = cp.bmat([[T_R, -T_I], [T_I, T_R]])
        constraints.append(big >> 0)

        self.constraints = constraints
        self.problem = cp.Problem(cp.Minimize(self.Omega), constraints)

    def solve(self, solver: str = "MOSEK", verbose: bool = False, **kwargs) -> V40Result:
        val = self.problem.solve(solver=solver, verbose=verbose, **kwargs)
        return V40Result(
            status=self.problem.status,
            Omega=float(self.Omega.value) if self.Omega.value is not None else float("nan"),
            primal=float(val) if val is not None else float("nan"),
        )


if __name__ == "__main__":
    for N in (16, 32):
        for K in (16, 32):
            try:
                prob = AutocorrLowerBoundV40(N=N, K=K, n_tests=33)
                out = prob.solve(solver="MOSEK", verbose=False)
                print(f"N={N:2d}, K={K:2d}: Omega={out.Omega:.6f}")
            except Exception as e:
                print(f"N={N}, K={K}: failed: {e}")
