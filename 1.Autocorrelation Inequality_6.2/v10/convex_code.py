"""v10: extend v9 to general (non-even) f using Fejer test at multiple points.

Key observation
---------------

For general (not necessarily even) admissible f, the autocorrelation
F = f*f is in general NOT even (counter-example: f = delta_{1/4} gives
f*f concentrated at 1/2).  The Fejer kernel can still be applied:

  (F_K * F)(t*) = int F_K(t* - s) F(s) ds  <=  ||F_K||_1 * ||F||_inf
                                            =  Omega.

We test this at multiple t* in [-1/2, 1/2] and lift the bilinear right-
hand side to a moment matrix M on (1, a_1, ..., a_K, b_1, ..., b_K).

For f real (general), Fourier coefficients a_k = int f cos(2 pi k x) dx,
b_k = int f sin(2 pi k x) dx.  Then  F_hat(k) = f_hat(k)^2 = (a_k - i b_k)^2
                                              = (a_k^2 - b_k^2) - 2 i a_k b_k.
Re F_hat(k) = a_k^2 - b_k^2,  -Im F_hat(k) = 2 a_k b_k.

The Fejer-mean of F at t* equals
  (F_K * F)(t*) = sum_{k} F_K_hat(k) F_hat(k) e^{2 pi i k t*}
                = 1 + 2 sum_{k=1}^K (1 - k/(K+1)) [
                       (a_k^2 - b_k^2) cos(2 pi k t*)
                       + 2 a_k b_k sin(2 pi k t*) ].

Lifting (a_k, b_k) to a moment matrix M of size 2K+1 with
M[0,0] = 1, M[0, k] = a_k, M[0, K+k] = b_k, the rank-one truth gives
M[k,k] = a_k^2, M[K+k, K+k] = b_k^2, M[k, K+k] = a_k b_k.  In the
relaxation, M >> 0 with these first-row identifications, the diagonal
entries are >= the corresponding squares, but cross terms are free.

Constraint at each t*:
  Omega >= 1 + 2 sum_{k=1}^K (1-k/(K+1)) [
              (M[k,k] - M[K+k,K+k]) cos(2 pi k t*)
              + 2 M[k, K+k] sin(2 pi k t*)
           ].

This is rigorous: the rank-one M^* = vv^T from any admissible f gives
LHS = the true Fejer mean at t*, which is <= ||F||_inf = Omega.

The relaxation may give a smaller Omega because M can be non-rank-one,
but it remains a valid lower bound.  We add the linear bounds on a_k,
b_k from the f-side discretization p as well, to pin down the
relationship.
"""

from __future__ import annotations

from dataclasses import dataclass

import cvxpy as cp
import numpy as np


@dataclass
class V10Result:
    status: str
    Omega: float
    primal: float


class AutocorrLowerBoundV10:
    def __init__(self, N: int = 16, K: int = 16, n_tests: int = 65) -> None:
        if N < 1 or K < 1:
            raise ValueError("N and K must be positive integers")
        self.N = N
        self.K = K
        self.dim_p = 2 * N
        self.L_f = 1.0 / (4 * N)
        self.n_tests = n_tests

        self.Omega = cp.Variable(nonneg=True, name="Omega")
        self.p = cp.Variable(self.dim_p, nonneg=True, name="p")
        self.a = cp.Variable(K + 1, name="a")
        self.b = cp.Variable(K + 1, name="b")
        self.M = cp.Variable((2 * K + 1, 2 * K + 1), symmetric=True, name="M")

        constraints: list = []

        # mass + a_0, b_0
        constraints.append(cp.sum(self.p) == 1)
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

        # f-side linear bounds for a_k, b_k from p
        cos_min = np.zeros((K + 1, self.dim_p))
        cos_max = np.zeros((K + 1, self.dim_p))
        sin_min = np.zeros((K + 1, self.dim_p))
        sin_max = np.zeros((K + 1, self.dim_p))
        for k in range(1, K + 1):
            for j in range(self.dim_p):
                a_left = -0.25 + j * self.L_f
                a_right = a_left + self.L_f
                xs = np.linspace(a_left, a_right, 401)
                cv = np.cos(2.0 * np.pi * k * xs)
                sv = np.sin(2.0 * np.pi * k * xs)
                cos_min[k, j] = cv.min()
                cos_max[k, j] = cv.max()
                sin_min[k, j] = sv.min()
                sin_max[k, j] = sv.max()
        for k in range(1, K + 1):
            constraints.append(self.a[k] >= cos_min[k] @ self.p)
            constraints.append(self.a[k] <= cos_max[k] @ self.p)
            constraints.append(self.b[k] >= sin_min[k] @ self.p)
            constraints.append(self.b[k] <= sin_max[k] @ self.p)

        # Fejer-kernel test at multiple t*
        ts = np.linspace(-0.5, 0.5, n_tests)
        weights = np.array([1.0 - k / (K + 1) for k in range(1, K + 1)])
        for t_star in ts:
            cos_vals = np.cos(2 * np.pi * np.arange(1, K + 1) * t_star)
            sin_vals = np.sin(2 * np.pi * np.arange(1, K + 1) * t_star)
            terms = []
            for k in range(1, K + 1):
                A_k = self.M[k, k] - self.M[K + k, K + k]  # lifted a_k^2 - b_k^2
                B_k = 2.0 * self.M[k, K + k]                # lifted 2 a_k b_k
                terms.append(weights[k - 1] * (A_k * cos_vals[k - 1] + B_k * sin_vals[k - 1]))
            expr = 1.0 + 2.0 * cp.sum(cp.hstack(terms))
            constraints.append(self.Omega >= expr)

        self.constraints = constraints
        self.problem = cp.Problem(cp.Minimize(self.Omega), constraints)

    def solve(self, solver: str = "MOSEK", verbose: bool = False, **kwargs) -> V10Result:
        val = self.problem.solve(solver=solver, verbose=verbose, **kwargs)
        return V10Result(
            status=self.problem.status,
            Omega=float(self.Omega.value) if self.Omega.value is not None else float("nan"),
            primal=float(val) if val is not None else float("nan"),
        )


if __name__ == "__main__":
    for N in (8, 16):
        for K in (8, 16):
            try:
                prob = AutocorrLowerBoundV10(N=N, K=K, n_tests=33)
                out = prob.solve(solver="MOSEK", verbose=False)
                print(f"N={N:2d}, K={K:2d}: Omega={out.Omega:.6f}")
            except Exception as e:
                print(f"N={N}, K={K}: failed: {e}")
