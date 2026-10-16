"""v32: v30 (multi-point Fejer) + v11 shifted windows (W) constraints.

Both v30 and v11 use rigorous, valid constraints that may add
non-trivial information.  v11 captures the (p, Q)-side bilinear
structure via window integrals; v30 captures the F-side peak via
pointwise Fejer means.  Combined, they should give a tighter LP
bound than either alone.
"""

from __future__ import annotations

from dataclasses import dataclass

import cvxpy as cp
import numpy as np


@dataclass
class V32Result:
    status: str
    Omega: float
    primal: float


class AutocorrLowerBoundV32:
    def __init__(self, N: int = 16, K: int = 16, n_tests: int = 51) -> None:
        self.N = N
        self.K = K
        self.dim = 2 * N
        self.L = 1.0 / (4 * N)
        self.n_tests = n_tests

        self.Omega = cp.Variable(nonneg=True, name="Omega")
        self.p = cp.Variable(self.dim, nonneg=True, name="p")
        self.Q = cp.Variable((self.dim, self.dim), symmetric=True, name="Q")
        self.a = cp.Variable(K + 1, name="a")
        self.b = cp.Variable(K + 1, name="b")
        self.M = cp.Variable((2 * K + 1, 2 * K + 1), symmetric=True, name="M")
        self.v = cp.Variable(K, nonneg=True, name="v")

        constraints: list = []

        block = cp.bmat([
            [np.array([[1.0]]), cp.reshape(self.p, (1, self.dim), order="C")],
            [cp.reshape(self.p, (self.dim, 1), order="C"), self.Q],
        ])
        constraints.append(block >> 0)
        constraints.append(cp.sum(self.p) == 1)
        constraints.append(self.Q >= 0)
        constraints.append(cp.sum(self.Q, axis=1) == self.p)
        constraints.append(self.a[0] == 1)
        constraints.append(self.b[0] == 0)

        constraints.append(self.M >> 0)
        constraints.append(self.M[0, 0] == 1)
        for k in range(1, K + 1):
            constraints.append(self.M[0, k] == self.a[k])
            constraints.append(self.M[k, 0] == self.a[k])
            constraints.append(self.M[0, K + k] == self.b[k])
            constraints.append(self.M[K + k, 0] == self.b[k])

        cos_min = np.zeros((K + 1, self.dim))
        cos_max = np.zeros((K + 1, self.dim))
        sin_min = np.zeros((K + 1, self.dim))
        sin_max = np.zeros((K + 1, self.dim))
        for k in range(1, K + 1):
            for j in range(self.dim):
                a_l = -0.25 + j * self.L
                a_r = a_l + self.L
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

        # v14 autocorrelation Fejer
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

        # v11 (W) shifted windows on Q
        widths = list(np.linspace(2.0 * self.L, 0.5, 8))
        offsets = [0.0] + [c * self.L for c in (1, 2, 4, 8)] + [-c * self.L for c in (1, 2, 4, 8)]
        for c in offsets:
            for h in widths:
                a_w = c - h
                b_w = c + h
                if a_w < -0.5 or b_w > 0.5:
                    continue
                lo = int(np.ceil((a_w + 0.5) / self.L)) + 2
                hi = int(np.floor((b_w + 0.5) / self.L))
                inside = []
                for j in range(1, self.dim + 1):
                    for kk in range(1, self.dim + 1):
                        if lo <= j + kk <= hi:
                            inside.append(self.Q[j - 1, kk - 1])
                if inside:
                    constraints.append(cp.sum(cp.hstack(inside)) <= 2.0 * h * self.Omega)

        self.constraints = constraints
        self.problem = cp.Problem(cp.Minimize(self.Omega), constraints)

    def solve(self, solver: str = "MOSEK", verbose: bool = False, **kwargs) -> V32Result:
        val = self.problem.solve(solver=solver, verbose=verbose, **kwargs)
        return V32Result(
            status=self.problem.status,
            Omega=float(self.Omega.value) if self.Omega.value is not None else float("nan"),
            primal=float(val) if val is not None else float("nan"),
        )


if __name__ == "__main__":
    for N in (16, 32, 50):
        for K in (16, 32):
            try:
                prob = AutocorrLowerBoundV32(N=N, K=K, n_tests=33)
                out = prob.solve(solver="MOSEK", verbose=False)
                print(f"N={N:2d}, K={K:2d}: Omega={out.Omega:.6f}")
            except Exception as e:
                print(f"N={N}, K={K}: failed: {e}")
