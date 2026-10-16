"""v17: parameter sweep of v14 at large N, K (per v14 verdict feedback #1).

Per the v9 verdict §5 recommendation, the right Fejer-style construct
for non-even f uses the *autocorrelation* F * F_check (where
F_check(t) = F(-t)) instead of F itself.  Key facts:

(A) (F * F_check)_hat(k) = F_hat(k) F_hat(-k) = F_hat(k) F_hat(k)*
                         = |F_hat(k)|^2 >= 0  (real, nonneg).
(B) (F * F_check) is also nonneg and has integral
    int (F * F_check) dt = (int F dt)^2 = 1.
(C) ||F * F_check||_inf = (F * F_check)(0) = int F^2 dt  (peak at 0).

The Fejer test:
   int F_K(t) (F * F_check)(t) dt  <=  ||F_K||_1 * ||F * F_check||_inf
                                    =  1 * int F^2.

Now use int F^2 <= ||F||_inf * int F = Omega:

   int F_K(t) (F * F_check)(t) dt  <=  Omega.

The LHS, by Parseval, equals
   sum_{|k|<=K} F_K_hat(k) * |F_hat(k)|^2
   = sum_{|k|<=K} (1 - |k|/(K+1)) * |F_hat(k)|^2
   = 1 + 2 sum_{k=1}^K (1 - k/(K+1)) * |F_hat(k)|^2

(using |F_hat(k)|^2 = |F_hat(-k)|^2 and F_hat(0) = 1).

The constraint becomes
   1 + 2 sum_{k=1}^K (1 - k/(K+1)) * |F_hat(k)|^2  <=  Omega.       (autoFejer)

The CRUCIAL feature: the coefficients of |F_hat(k)|^2 are POSITIVE
(Fejer weights are nonneg).  In the SDP lift, we use
|F_hat(k)|^2 = |f_hat(k)|^4 = (a_k^2 + b_k^2)^2, and we lift this
via two layers:

(L1) Introduce u_k >= 0 with the LIFTED PSD bound
     u_k <= M[k,k] + M[K+k, K+k]                         -- wrong direction
  Actually we want u_k = a_k^2 + b_k^2 = |f_hat(k)|^2 in the truth.
  In the lift, u_k >= a_k^2 + b_k^2 isn't directly available, since
  M[k,k] >= a_k^2 only gives a LOWER bound on M[k,k] (not on a_k^2).

  Take u_k := M[k,k] + M[K+k, K+k] (linear in M).
  Then in the rank-one truth, u_k = a_k^2 + b_k^2.  In the relaxation,
  u_k could be > a_k^2 + b_k^2 (since M's diagonal can be inflated).
  This is the OPPOSITE direction of what we want for a tight bound.

  But for our constraint, the relaxation will set M[k,k] and M[K+k,K+k]
  to their MINIMUM values to minimize u_k (since the constraint is
  sum (positive weight) * u_k^2 <= Omega).  At the minimum,
  u_k = a_k^2 + b_k^2 (the rank-one truth value), as long as the
  PSD lift is tight.

(L2) Introduce v_k >= 0 with v_k >= u_k^2 (rotated SOC, convex).
     In the rank-one truth, v_k = u_k^2 = (a_k^2 + b_k^2)^2 = |F_hat(k)|^2.

Constraint (autoFejer) becomes
   Omega  >=  1 + 2 sum_{k=1}^K (1 - k/(K+1)) * v_k.

Validity: at the rank-one truth M = vv^T, M[k,k] = a_k^2,
M[K+k, K+k] = b_k^2, so u_k = a_k^2 + b_k^2 (= |f_hat(k)|^2 truth).
v_k can be set to u_k^2 (the smallest value satisfying the SOC),
giving v_k = (a_k^2+b_k^2)^2 = |F_hat(k)|^2.  The constraint then
matches (autoFejer) at equality.

The relaxation may have u_k > truth and v_k > u_k^2, but the LP
minimizes Omega, so it sets v_k = u_k^2 (the tightest), and u_k =
M[k,k] + M[K+k, K+k] >= a_k^2 + b_k^2 = truth (the tightest given
PSD).  So the LP picks the right values, and the bound on Omega is
the rank-one truth bound, which is non-trivial.

Numerical performance is reported in rigorousproof.md.
"""

from __future__ import annotations

from dataclasses import dataclass

import cvxpy as cp
import numpy as np


@dataclass
class V14Result:
    status: str
    Omega: float
    primal: float


class AutocorrLowerBoundV14:
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
        self.v = cp.Variable(K, nonneg=True, name="v")  # |F_hat(k)|^2 lifted, k=1..K

        constraints: list = []

        # mass and basics
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

        # Linear bounds on a_k, b_k from p
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

        # SOC: v_k >= (M[k,k] + M[K+k, K+k])^2 (rotated SOC, convex).
        # Equivalent: v_k * 1 >= (M[k,k] + M[K+k, K+k])^2.
        # Use cp.quad_over_lin: cp.quad_over_lin(M[k,k]+M[K+k,K+k], 1) <= v_k.
        for k in range(1, K + 1):
            sum_diag = self.M[k, k] + self.M[K + k, K + k]
            constraints.append(self.v[k - 1] >= cp.square(sum_diag))

        # Fejer autocorrelation constraint:
        # Omega >= 1 + 2 * sum_{k=1}^K (1 - k/(K+1)) * v_{k-1}
        weights = np.array([1.0 - k / (K + 1) for k in range(1, K + 1)])
        constraints.append(self.Omega >= 1.0 + 2.0 * (weights @ self.v))

        self.constraints = constraints
        self.problem = cp.Problem(cp.Minimize(self.Omega), constraints)

    def solve(self, solver: str = "MOSEK", verbose: bool = False, **kwargs) -> V14Result:
        val = self.problem.solve(solver=solver, verbose=verbose, **kwargs)
        return V14Result(
            status=self.problem.status,
            Omega=float(self.Omega.value) if self.Omega.value is not None else float("nan"),
            primal=float(val) if val is not None else float("nan"),
        )


if __name__ == "__main__":
    import time
    for N, K in [(50, 32), (50, 64), (75, 64), (100, 64), (100, 128), (150, 128)]:
        t = time.time()
        try:
            prob = AutocorrLowerBoundV14(N=N, K=K)
            out = prob.solve(solver="MOSEK", verbose=False)
            print(f"N={N:3d}, K={K:3d}: Omega={out.Omega:.6f}  ({time.time()-t:.1f}s)")
        except Exception as e:
            print(f"N={N}, K={K}: failed: {e}")
