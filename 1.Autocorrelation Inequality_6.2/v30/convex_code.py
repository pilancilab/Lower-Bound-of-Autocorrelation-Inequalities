"""v30: pointwise Fejer mean tests at multiple t*.

Per v29's diagnostic: v14 captures int F^2 (the average), not max F
(the peak).  To get closer to max F, use the **Fejer mean of F at
multiple points t***:
   sigma_K(t*) := (F * F_K)(t*),
which approximates F(t*) for large K (with error decaying as 1/K).

The constraint:
   sigma_K(t*) <= ||F||_inf = Omega.

The lifted form involves both Re F_hat(k) and Im F_hat(k):
   sigma_K(t*) = sum_{|k|<=K} F_K_hat(k) * F_hat(k) e^{2 pi i k t*}
              = 1 + 2 sum_{k=1}^{K} (1 - k/(K+1)) * [Re F_hat(k) cos(2 pi k t*)
                                                     + Im F_hat(k) (-sin(2 pi k t*))]
              = 1 + 2 sum_k (1 - k/(K+1)) * [(M[k,k] - M[K+k, K+k]) cos(2 pi k t*)
                                             + 2 M[k, K+k] sin(2 pi k t*)]
                                             # since Im F_hat = -2 a b => -Im... cos
                                             # let me re-derive the sign

For real F: F_hat(-k) = conj F_hat(k). The Fourier inversion at t*:
   F(t*) = sum_k F_hat(k) e^{2 pi i k t*}
         = F_hat(0) + sum_{k>=1} [F_hat(k) e^{2 pi i k t*} + F_hat(-k) e^{-2 pi i k t*}]
         = 1 + 2 sum_{k>=1} Re[F_hat(k) e^{2 pi i k t*}]
         = 1 + 2 sum_{k>=1} [Re F_hat(k) cos(2 pi k t*) - Im F_hat(k) sin(2 pi k t*)].

Sigma_K(t*) is the same with weights (1 - k/(K+1)).

In our lift: Re F_hat(k) = M[k,k] - M[K+k, K+k], Im F_hat(k) = -2 M[k, K+k].

So sigma_K(t*) = 1 + 2 sum_{k=1}^K (1 - k/(K+1)) * [
                    (M[k,k] - M[K+k, K+k]) cos(2 pi k t*)
                    + 2 M[k, K+k] sin(2 pi k t*)  # since -Im = +2 a b
                  ].

Constraint:  Omega >= sigma_K(t*).

Note: this is the SIGNED (not autocorrelation) version, and is
SUBJECT TO the same diagonal-inflation issue as v2.  However, by
COMBINING with the v14 autocorrelation Fejer (which is unsigned),
we may get a useful bound.
"""

from __future__ import annotations

from dataclasses import dataclass

import cvxpy as cp
import numpy as np


@dataclass
class V30Result:
    status: str
    Omega: float
    primal: float


class AutocorrLowerBoundV30:
    def __init__(self, N: int = 8, K: int = 8, n_tests: int = 33) -> None:
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

        # v14 autocorrelation Fejer
        for k in range(1, K + 1):
            sum_diag = self.M[k, k] + self.M[K + k, K + k]
            constraints.append(self.v[k - 1] >= cp.square(sum_diag))
        weights = np.array([1.0 - k / (K + 1) for k in range(1, K + 1)])
        constraints.append(self.Omega >= 1.0 + 2.0 * (weights @ self.v))

        # NEW v30: pointwise Fejer mean tests at multiple t*
        ts = np.linspace(-0.5, 0.5, n_tests)
        for t_star in ts:
            cos_vals = np.cos(2 * np.pi * np.arange(1, K + 1) * t_star)
            sin_vals = np.sin(2 * np.pi * np.arange(1, K + 1) * t_star)
            terms = []
            for k in range(1, K + 1):
                # Re F_hat = M[k,k] - M[K+k, K+k]
                # -Im F_hat = 2 M[k, K+k]
                re_F = self.M[k, k] - self.M[K + k, K + k]
                neg_im_F = 2.0 * self.M[k, K + k]  # this is -Im F_hat
                weight_k = 1.0 - k / (K + 1)
                terms.append(weight_k * (re_F * cos_vals[k - 1] + neg_im_F * sin_vals[k - 1]))
            sigma_K_t = 1.0 + 2.0 * cp.sum(cp.hstack(terms))
            constraints.append(self.Omega >= sigma_K_t)

        self.constraints = constraints
        self.problem = cp.Problem(cp.Minimize(self.Omega), constraints)

    def solve(self, solver: str = "MOSEK", verbose: bool = False, **kwargs) -> V30Result:
        val = self.problem.solve(solver=solver, verbose=verbose, **kwargs)
        return V30Result(
            status=self.problem.status,
            Omega=float(self.Omega.value) if self.Omega.value is not None else float("nan"),
            primal=float(val) if val is not None else float("nan"),
        )


if __name__ == "__main__":
    for N in (8, 16):
        for K in (8, 16):
            try:
                prob = AutocorrLowerBoundV30(N=N, K=K, n_tests=33)
                out = prob.solve(solver="MOSEK", verbose=False)
                print(f"N={N:2d}, K={K:2d}: Omega={out.Omega:.6f}")
            except Exception as e:
                print(f"N={N}, K={K}: failed: {e}")
