"""v24: v14 + f-side L^infty variable M.

Introduce LP variable M >= 0 representing an upper bound on the
density of f.  The constraints encode that M is at least the
discrete cell density p_j / L_f for each j, AND M is at least Omega
(via int f^2 <= ||f||_inf * int f = M and int f^2 >= ||f*f||_inf =
Omega from Cauchy-Schwarz).

The new constraints:
* M_inf >= p_j / L_f  for all j  (f's L_inf >= cell average)
* M_inf >= Omega       (from int f^2 >= Omega and int f^2 <= M_inf)

These are LINEAR.  The first encodes that the LP can't have a
single cell with mass exceeding M_inf * L_f.  The second couples
M_inf to Omega.

Note: M_inf >= Omega is a CONSEQUENCE of int f^2 >= Omega and
int f^2 <= M_inf, both rigorous.

This is rigorously valid: every admissible f gives M_inf = ||f||_inf
which satisfies both bounds.
"""

from __future__ import annotations

from dataclasses import dataclass

import cvxpy as cp
import numpy as np


@dataclass
class V24Result:
    status: str
    Omega: float
    primal: float


class AutocorrLowerBoundV24:
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
        self.M_inf = cp.Variable(nonneg=True, name="M_inf")

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

        # NEW v24: M_inf bounds
        for j in range(self.dim_p):
            constraints.append(self.M_inf >= self.p[j] / self.L_f)
        constraints.append(self.M_inf >= self.Omega)
        # Also: M_inf >= 2 (from ||f||_inf >= 2 derived from Cauchy-Schwarz on f)
        # Proof: int f^4 >= (int f^2)^2 / (1/2) = 2 (int f^2)^2 >= 2 * 4 = 8.
        # And int f^4 <= ||f||_inf^2 * int f^2 <= M_inf^3 (using int f^2 <= M_inf).
        # So 8 <= M_inf^3 -> M_inf >= 2.
        constraints.append(self.M_inf >= 2.0)

        self.constraints = constraints
        self.problem = cp.Problem(cp.Minimize(self.Omega), constraints)

    def solve(self, solver: str = "MOSEK", verbose: bool = False, **kwargs) -> V24Result:
        val = self.problem.solve(solver=solver, verbose=verbose, **kwargs)
        return V24Result(
            status=self.problem.status,
            Omega=float(self.Omega.value) if self.Omega.value is not None else float("nan"),
            primal=float(val) if val is not None else float("nan"),
        )


if __name__ == "__main__":
    for N in (8, 16):
        for K in (8, 16):
            try:
                prob = AutocorrLowerBoundV24(N=N, K=K)
                out = prob.solve(solver="MOSEK", verbose=False)
                print(f"N={N:2d}, K={K:2d}: Omega={out.Omega:.6f}")
            except Exception as e:
                print(f"N={N}, K={K}: failed: {e}")
