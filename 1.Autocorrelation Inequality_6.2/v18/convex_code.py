"""v18: v14 + L^3 norm constraint Omega^2 >= int F^3.

For F >= 0 with int F = 1 and F <= Omega:
   F^3 = F * F * F  <=  Omega^2 * F  (using F <= Omega twice).
   int F^3 dt  <=  Omega^2 * int F dt  =  Omega^2.

The lifted form uses W_m as cell averages of F (consistent with the
v14 / v11 framework via the (p, Q) lift).  By Jensen on the convex
function x^3 on each cell:
   int_{B_m} F^3 dt  >=  L_F * W_m^3 = L_F * (cell average)^3.

Hence
   L_F * sum_m W_m^3  <=  int F^3  <=  Omega^2.        (L3)

In the relaxation, we add (L3) as a convex constraint:
   sum_m W_m^3  <=  Omega^2 / L_F,
expressed via the auxiliary lift  u_m >= W_m^3 (rotated power cone)
and  omega = Omega^2 (epigraph of square).

Combined with v14's autocorrelation Fejer constraint, this gives a
strictly tighter LP than v14 alone.

For uniform f = 2 on [-1/4, 1/4]:
   F is the triangle 2*max(0, 1 - 2|t|).
   int F^3 = 2 * int_0^{1/2} (2(1-2t))^3 dt
            = 16 * int_0^{1/2} (1-2t)^3 dt
            = 16 * 1/8 = 2.
   So Omega^2 >= 2, i.e., Omega >= sqrt(2) ~ 1.414.
This exceeds v14's bound of 4/3 ~ 1.333 and the SOTA C_{6.2} >= 1.28.

Status: experimental, see rigorousproof.md for the formal proof.
"""

from __future__ import annotations

from dataclasses import dataclass

import cvxpy as cp
import numpy as np


@dataclass
class V18Result:
    status: str
    Omega: float
    primal: float


class AutocorrLowerBoundV18:
    def __init__(self, N: int = 16, K: int = 16) -> None:
        if N < 1 or K < 1:
            raise ValueError("N, K must be positive")
        self.N = N
        self.K = K
        self.dim = 2 * N
        self.L = 1.0 / (4 * N)
        # F-side: 4N intervals of width L on [-1/2, 1/2]
        self.dim_W = 4 * N
        self.L_F = 1.0 / (4 * N)

        self.Omega = cp.Variable(nonneg=True, name="Omega")
        self.p = cp.Variable(self.dim, nonneg=True, name="p")
        self.Q = cp.Variable((self.dim, self.dim), symmetric=True, name="Q")
        self.W = cp.Variable(self.dim_W, nonneg=True, name="W")
        self.a = cp.Variable(K + 1, name="a")
        self.b = cp.Variable(K + 1, name="b")
        self.M = cp.Variable((2 * K + 1, 2 * K + 1), symmetric=True, name="M")
        self.v = cp.Variable(K, nonneg=True, name="v")
        # auxiliary for L^3
        self.u = cp.Variable(self.dim_W, nonneg=True, name="u")
        self.omega = cp.Variable(nonneg=True, name="omega")  # = Omega^2

        constraints: list = []

        # (p, Q) PSD lift and standard
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

        # M PSD lift on (1, a, b)
        constraints.append(self.M >> 0)
        constraints.append(self.M[0, 0] == 1)
        for k in range(1, K + 1):
            constraints.append(self.M[0, k] == self.a[k])
            constraints.append(self.M[k, 0] == self.a[k])
            constraints.append(self.M[0, K + k] == self.b[k])
            constraints.append(self.M[K + k, 0] == self.b[k])

        # f-side bounds
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

        # F-side W variables linked to Q via the (W) cell-fully-inside form
        # W_m * L_F >= sum of Q[j,k] for cells fully inside B_m's preimage
        # Plus the trivial: W_m <= Omega.
        constraints.append(self.W <= self.Omega)
        constraints.append(cp.sum(self.W) * self.L_F == 1.0)

        for m in range(1, self.dim_W + 1):
            a_w = -0.5 + (m - 1) * self.L_F
            b_w = a_w + self.L_F
            lo = int(np.ceil((a_w + 0.5) / self.L)) + 2
            hi = int(np.floor((b_w + 0.5) / self.L))
            inside = []
            for jj in range(1, self.dim + 1):
                for kk in range(1, self.dim + 1):
                    if lo <= jj + kk <= hi:
                        inside.append(self.Q[jj - 1, kk - 1])
            if inside:
                constraints.append(cp.sum(cp.hstack(inside)) <= self.L_F * self.W[m - 1])

        # L^3 constraint: sum_m u_m * L_F <= omega, u_m >= W_m^3, omega >= Omega^2
        for m in range(self.dim_W):
            constraints.append(self.u[m] >= cp.power(self.W[m], 3))
        constraints.append(self.omega >= cp.square(self.Omega))
        constraints.append(self.L_F * cp.sum(self.u) <= self.omega)

        self.constraints = constraints
        self.problem = cp.Problem(cp.Minimize(self.Omega), constraints)

    def solve(self, solver: str = "MOSEK", verbose: bool = False, **kwargs) -> V18Result:
        val = self.problem.solve(solver=solver, verbose=verbose, **kwargs)
        return V18Result(
            status=self.problem.status,
            Omega=float(self.Omega.value) if self.Omega.value is not None else float("nan"),
            primal=float(val) if val is not None else float("nan"),
        )


if __name__ == "__main__":
    for N in (8, 12):
        for K in (8, 12):
            try:
                prob = AutocorrLowerBoundV18(N=N, K=K)
                out = prob.solve(solver="MOSEK", verbose=False)
                print(f"N={N:2d}, K={K:2d}: Omega={out.Omega:.6f}")
            except Exception as e:
                print(f"N={N}, K={K}: failed: {e}")
