"""v22: v14 with tighter midpoint quadrature for a_k, b_k cell bounds.

Per v20 verdict's Comment 3: replace the loose min/max cell-bound on
a_k = int f cos(2 pi k x) dx with a midpoint quadrature plus explicit
second-order error.

Specifically, for f >= 0 with int_{I_j} f = p_j, the midpoint rule
on I_j gives
   int_{I_j} f * cos(2 pi k x) dx
       =  p_j * cos(2 pi k * x_j)  +  R_{j, k}
where x_j is the midpoint of I_j, and the remainder R_{j, k} is
bounded by

   |R_{j, k}|  <=  (L_f^2 / 8) * (2 pi k)^2 * p_j  =  (pi^2 k^2 L_f^2 / 2) * p_j.

Sum over j:
   |a_k - sum_j p_j cos(2 pi k x_j)|  <=  (pi^2 k^2 L_f^2 / 2).

This is *much* tighter than the loose bound when k is small relative
to 1/L_f.

For uniform p (= 1/dim_p), the midpoint sum approaches
   int_{-1/4}^{1/4} 2 cos(2 pi k x) dx  =  2 sin(pi k / 2) / (pi k),
which is the true a_k for f = 2 (uniform).  Good agreement.

The lifted constraint becomes:
   sum_j p_j cos(2 pi k x_j)  -  err_k  <=  a_k  <=  sum_j p_j cos(2 pi k x_j)  +  err_k,
where err_k = pi^2 k^2 L_f^2 / 2.
"""

from __future__ import annotations

from dataclasses import dataclass

import cvxpy as cp
import numpy as np


@dataclass
class V22Result:
    status: str
    Omega: float
    primal: float


class AutocorrLowerBoundV22:
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

        # midpoint quadrature for a_k, b_k bounds with explicit error
        # x_j = -1/4 + (j - 0.5) * L_f
        x_mid = np.array([-0.25 + (j - 0.5) * self.L_f for j in range(1, self.dim_p + 1)])
        cos_mid = np.zeros((K + 1, self.dim_p))
        sin_mid = np.zeros((K + 1, self.dim_p))
        for k in range(1, K + 1):
            cos_mid[k] = np.cos(2 * np.pi * k * x_mid)
            sin_mid[k] = np.sin(2 * np.pi * k * x_mid)
        for k in range(1, K + 1):
            err_k = (np.pi ** 2) * (k ** 2) * (self.L_f ** 2) / 2.0
            constraints.append(self.a[k] >= cos_mid[k] @ self.p - err_k)
            constraints.append(self.a[k] <= cos_mid[k] @ self.p + err_k)
            constraints.append(self.b[k] >= sin_mid[k] @ self.p - err_k)
            constraints.append(self.b[k] <= sin_mid[k] @ self.p + err_k)

        for k in range(1, K + 1):
            sum_diag = self.M[k, k] + self.M[K + k, K + k]
            constraints.append(self.v[k - 1] >= cp.square(sum_diag))

        weights = np.array([1.0 - k / (K + 1) for k in range(1, K + 1)])
        constraints.append(self.Omega >= 1.0 + 2.0 * (weights @ self.v))

        self.constraints = constraints
        self.problem = cp.Problem(cp.Minimize(self.Omega), constraints)

    def solve(self, solver: str = "MOSEK", verbose: bool = False, **kwargs) -> V22Result:
        val = self.problem.solve(solver=solver, verbose=verbose, **kwargs)
        return V22Result(
            status=self.problem.status,
            Omega=float(self.Omega.value) if self.Omega.value is not None else float("nan"),
            primal=float(val) if val is not None else float("nan"),
        )


if __name__ == "__main__":
    for N in (8, 16, 32):
        for K in (8, 16, 32):
            try:
                prob = AutocorrLowerBoundV22(N=N, K=K)
                out = prob.solve(solver="MOSEK", verbose=False)
                print(f"N={N:2d}, K={K:2d}: Omega={out.Omega:.6f}")
            except Exception as e:
                print(f"N={N}, K={K}: failed: {e}")
