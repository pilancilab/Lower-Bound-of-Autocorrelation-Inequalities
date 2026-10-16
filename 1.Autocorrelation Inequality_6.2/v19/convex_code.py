"""v19: v14 with autocorrelation Fejer constraints at multiple K values.

The v14 autocorrelation Fejer constraint
   Omega >= 1 + 2 sum_{k=1}^K (1 - k/(K+1)) v_k
holds for any K.  Different K give different (linearly independent)
constraints, all valid.  v19 adds the constraint at multiple K values
simultaneously.
"""

from __future__ import annotations

from dataclasses import dataclass

import cvxpy as cp
import numpy as np


@dataclass
class V19Result:
    status: str
    Omega: float
    primal: float


class AutocorrLowerBoundV19:
    def __init__(self, N: int = 8, K_max: int = 16) -> None:
        if N < 1 or K_max < 1:
            raise ValueError("N, K_max must be positive")
        self.N = N
        self.K_max = K_max
        self.dim_p = 2 * N
        self.L_f = 1.0 / (4 * N)

        self.Omega = cp.Variable(nonneg=True, name="Omega")
        self.p = cp.Variable(self.dim_p, nonneg=True, name="p")
        self.a = cp.Variable(K_max + 1, name="a")
        self.b = cp.Variable(K_max + 1, name="b")
        self.M = cp.Variable((2 * K_max + 1, 2 * K_max + 1), symmetric=True, name="M")
        self.v = cp.Variable(K_max, nonneg=True, name="v")

        constraints: list = []

        constraints.append(cp.sum(self.p) == 1)
        constraints.append(self.a[0] == 1)
        constraints.append(self.b[0] == 0)

        constraints.append(self.M >> 0)
        constraints.append(self.M[0, 0] == 1)
        for k in range(1, K_max + 1):
            constraints.append(self.M[0, k] == self.a[k])
            constraints.append(self.M[k, 0] == self.a[k])
            constraints.append(self.M[0, K_max + k] == self.b[k])
            constraints.append(self.M[K_max + k, 0] == self.b[k])

        cos_min = np.zeros((K_max + 1, self.dim_p))
        cos_max = np.zeros((K_max + 1, self.dim_p))
        sin_min = np.zeros((K_max + 1, self.dim_p))
        sin_max = np.zeros((K_max + 1, self.dim_p))
        for k in range(1, K_max + 1):
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
        for k in range(1, K_max + 1):
            constraints.append(self.a[k] >= cos_min[k] @ self.p)
            constraints.append(self.a[k] <= cos_max[k] @ self.p)
            constraints.append(self.b[k] >= sin_min[k] @ self.p)
            constraints.append(self.b[k] <= sin_max[k] @ self.p)

        for k in range(1, K_max + 1):
            sum_diag = self.M[k, k] + self.M[K_max + k, K_max + k]
            constraints.append(self.v[k - 1] >= cp.square(sum_diag))

        # Fejer constraint at MULTIPLE K values
        for K_eff in range(2, K_max + 1, 2):
            weights = np.array([1.0 - k / (K_eff + 1) for k in range(1, K_eff + 1)])
            v_slice = cp.hstack([self.v[k - 1] for k in range(1, K_eff + 1)])
            constraints.append(self.Omega >= 1.0 + 2.0 * (weights @ v_slice))

        self.constraints = constraints
        self.problem = cp.Problem(cp.Minimize(self.Omega), constraints)

    def solve(self, solver: str = "MOSEK", verbose: bool = False, **kwargs) -> V19Result:
        val = self.problem.solve(solver=solver, verbose=verbose, **kwargs)
        return V19Result(
            status=self.problem.status,
            Omega=float(self.Omega.value) if self.Omega.value is not None else float("nan"),
            primal=float(val) if val is not None else float("nan"),
        )


if __name__ == "__main__":
    for N in (8, 16, 32):
        for Kmax in (16, 32):
            try:
                prob = AutocorrLowerBoundV19(N=N, K_max=Kmax)
                out = prob.solve(solver="MOSEK", verbose=False)
                print(f"N={N:2d}, Kmax={Kmax:2d}: Omega={out.Omega:.6f}")
            except Exception as e:
                print(f"N={N}, Kmax={Kmax}: failed: {e}")
