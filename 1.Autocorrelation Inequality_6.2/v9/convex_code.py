"""v9: rigorous Fejer-kernel test combined with v7 (W) constraints.

Key insight from v5 verdict's Comment 2 (suggested test functions):
For any nonneg trig polynomial P(t) on [-1/2, 1/2],
   integral P(t) F(t) dt  <=  Omega * integral P(t) dt.
This is rigorous (no Gibbs/tail issue) provided P >= 0 pointwise.

We use the Fejer kernel of order K:
   F_K(t) = sum_{|k|<=K} (1 - |k|/(K+1)) * exp(2 pi i k t),
which is nonneg with int F_K dt = 1.  Then
   integral F_K(t) F(t) dt = sum_{k} (1 - |k|/(K+1)) * F_hat(k)
                            = 1 + 2 sum_{k=1}^K (1 - k/(K+1)) * Re F_hat(k).
For F = f*f and even f, Re F_hat(k) = a_k^2 (with b_k = 0).

Hence the rigorous constraint
   1 + 2 sum_{k=1}^K (1 - k/(K+1)) * a_k^2  <=  Omega * 1.

Lifting a_k^2 to M[k,k] >= a_k^2 (PSD lift) gives the LP constraint
   Omega  >=  1 + 2 sum_{k=1}^K (1 - k/(K+1)) * M[k,k]
            (because the lifted M[k,k] >= a_k^2 makes the lifted RHS
             *larger* than the true RHS, so the inequality is preserved
             in the wrong direction... wait this needs more thought.)

Actually, the SAFE direction:
   Omega  >=  true RHS  =  1 + 2 sum_k (1 - k/(K+1)) a_k^2.
The lifted M[k,k] is >= a_k^2, so we have
   1 + 2 sum_k (1 - k/(K+1)) M[k,k]  >=  1 + 2 sum_k (1 - k/(K+1)) a_k^2  =  true RHS.
So enforcing  Omega >= 1 + 2 sum (1-k/(K+1)) M[k,k]  enforces
   Omega >= true RHS, which is the desired rigorous inequality.

So the constraint  Omega >= 1 + 2 sum (1 - k/(K+1)) M[k,k]  is rigorously
valid for every admissible even f with the rank-one lift M = (1, a_1..a_K)
* (1, a_1..a_K)^T.

Note this only bounds C^{even}_{6.2}.

The relaxation must also enforce M >> 0 with M[0,0]=1 and M[0,k]=a_k,
plus the linear bounds on a_k from p (the f-side discretization).
"""

from __future__ import annotations

from dataclasses import dataclass

import cvxpy as cp
import numpy as np


@dataclass
class V9Result:
    status: str
    Omega: float
    primal: float


class AutocorrLowerBoundV9:
    def __init__(self, N: int = 8, K: int = 12) -> None:
        if N < 1 or K < 1:
            raise ValueError("N and K must be positive integers")
        self.N = N
        self.K = K
        self.dim_p = 2 * N
        self.L_f = 1.0 / (4 * N)

        self.Omega = cp.Variable(nonneg=True, name="Omega")
        self.p = cp.Variable(self.dim_p, nonneg=True, name="p")
        self.a = cp.Variable(K + 1, name="a")
        self.M = cp.Variable((K + 1, K + 1), symmetric=True, name="M")

        constraints: list = []

        # mass + a_0
        constraints.append(cp.sum(self.p) == 1)
        constraints.append(self.a[0] == 1)

        # even-f symmetry
        for i in range(self.dim_p // 2):
            constraints.append(self.p[i] == self.p[self.dim_p - 1 - i])

        # PSD lift on M
        constraints.append(self.M >> 0)
        constraints.append(self.M[0, 0] == 1)
        for k in range(1, K + 1):
            constraints.append(self.M[0, k] == self.a[k])
            constraints.append(self.M[k, 0] == self.a[k])

        # Linear bounds on a_k from p (f-side; cells of width L_f on [-1/4, 1/4])
        cos_min = np.zeros((K + 1, self.dim_p))
        cos_max = np.zeros((K + 1, self.dim_p))
        for k in range(1, K + 1):
            for j in range(self.dim_p):
                a_left = -0.25 + j * self.L_f
                a_right = a_left + self.L_f
                xs = np.linspace(a_left, a_right, 401)
                cv = np.cos(2.0 * np.pi * k * xs)
                cos_min[k, j] = cv.min()
                cos_max[k, j] = cv.max()
        for k in range(1, K + 1):
            constraints.append(self.a[k] >= cos_min[k] @ self.p)
            constraints.append(self.a[k] <= cos_max[k] @ self.p)

        # Fejer kernel constraint
        # Omega >= 1 + 2 sum_{k=1}^K (1 - k/(K+1)) M[k,k]
        weights = np.array([1.0 - k / (K + 1) for k in range(1, K + 1)])
        m_diags = cp.hstack([self.M[k, k] for k in range(1, K + 1)])
        constraints.append(self.Omega >= 1.0 + 2.0 * (weights @ m_diags))

        self.constraints = constraints
        self.problem = cp.Problem(cp.Minimize(self.Omega), constraints)

    def solve(self, solver: str = "MOSEK", verbose: bool = False, **kwargs) -> V9Result:
        val = self.problem.solve(solver=solver, verbose=verbose, **kwargs)
        return V9Result(
            status=self.problem.status,
            Omega=float(self.Omega.value) if self.Omega.value is not None else float("nan"),
            primal=float(val) if val is not None else float("nan"),
        )


if __name__ == "__main__":
    for N in (8, 16):
        for K in (4, 8, 16, 32):
            try:
                prob = AutocorrLowerBoundV9(N=N, K=K)
                out = prob.solve(solver="MOSEK", verbose=False)
                print(f"N={N:2d}, K={K:2d}: Omega={out.Omega:.6f}")
            except Exception as exc:
                print(f"N={N}, K={K}: failed: {exc}")
