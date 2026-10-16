"""v20: v14 with Lasserre level-2 lift on (a, b).

Per v15 verdict: strengthen the SDP lift on M from degree 1 to
degree 2 in the variables (a_1, ..., a_K, b_1, ..., b_K).  The
level-2 moment matrix has rows/cols indexed by all monomials of
degree <= 2 in these 2K variables.

The structure is:
* y_alpha for multi-indices alpha = (alpha_a, alpha_b) with
  alpha_a in N^K, alpha_b in N^K, |alpha| := |alpha_a| + |alpha_b| <= 4.
* M2[i, j] = y_{alpha_i + alpha_j}, where alpha_i ranges over
  multi-indices of degree <= 2.
* PSD constraint on M2.

For K = 4 (8 variables), # multi-indices of degree <= 2 is
C(10, 2) = 45.  M2 is 45x45.

We then identify:
* y_0 = 1
* y_{e_k^a} = a_k                               (linear)
* y_{e_k^b} = b_k                               (linear)
* y_{e_k^a + e_l^a} = a_k a_l                   (lifted)
* y_{e_k^b + e_l^b} = b_k b_l                   (lifted)
* y_{e_k^a + e_l^b} = a_k b_l                   (lifted)

The v14 lifted variables M[k,k] = a_k^2 = y_{2 e_k^a},
M[K+k, K+k] = b_k^2 = y_{2 e_k^b},
M[k, K+k] = a_k b_k = y_{e_k^a + e_k^b}.

The Fejer constraint and rotated SOC for v_k are unchanged.

The PSD constraint on M2 (which has more rows/cols than the v14 M)
gives strictly more constraints, hence a tighter relaxation.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations_with_replacement

import cvxpy as cp
import numpy as np


@dataclass
class V20Result:
    status: str
    Omega: float
    primal: float


def enum_alphas(n: int, max_deg: int):
    """Enumerate all multi-indices of degree <= max_deg in n variables."""
    out = []
    for d in range(max_deg + 1):
        for combo in combinations_with_replacement(range(n), d):
            alpha = [0] * n
            for i in combo:
                alpha[i] += 1
        # Note: combinations_with_replacement gives sorted tuples; convert to multi-index.
            out.append(tuple(alpha))
    return out


def enum_alphas_full(n: int, max_deg: int):
    out = []
    for d in range(max_deg + 1):
        for combo in combinations_with_replacement(range(n), d):
            alpha = [0] * n
            for i in combo:
                alpha[i] += 1
            out.append(tuple(alpha))
    return out


class AutocorrLowerBoundV20:
    def __init__(self, N: int = 8, K: int = 4) -> None:
        if N < 1 or K < 1:
            raise ValueError("N, K must be positive")
        self.N = N
        self.K = K
        self.dim_p = 2 * N
        self.L_f = 1.0 / (4 * N)

        # 2K variables: a_1..a_K, b_1..b_K, indexed 0..2K-1
        n_vars = 2 * K
        alphas = enum_alphas_full(n_vars, 4)
        self.alpha_to_idx = {a: i for i, a in enumerate(alphas)}
        self.alphas = alphas

        self.Omega = cp.Variable(nonneg=True, name="Omega")
        self.p = cp.Variable(self.dim_p, nonneg=True, name="p")
        self.a = cp.Variable(K + 1, name="a")
        self.b = cp.Variable(K + 1, name="b")
        self.y = cp.Variable(len(alphas), name="y")
        self.v = cp.Variable(K, nonneg=True, name="v")

        constraints: list = []

        # mass and a_0, b_0
        constraints.append(cp.sum(self.p) == 1)
        constraints.append(self.a[0] == 1)
        constraints.append(self.b[0] == 0)

        # y identification with monomials of degree 0 and 1
        zero_alpha = tuple([0] * n_vars)
        constraints.append(self.y[self.alpha_to_idx[zero_alpha]] == 1)
        for k in range(K):
            ek_a = [0] * n_vars
            ek_a[k] = 1
            constraints.append(self.y[self.alpha_to_idx[tuple(ek_a)]] == self.a[k + 1])
            ek_b = [0] * n_vars
            ek_b[K + k] = 1
            constraints.append(self.y[self.alpha_to_idx[tuple(ek_b)]] == self.b[k + 1])

        # Build moment matrix M2 of size {alphas with deg <= 2}
        alphas2 = [a for a in alphas if sum(a) <= 2]
        size2 = len(alphas2)
        rows = []
        for i in range(size2):
            row = []
            for j in range(size2):
                summed = tuple(alphas2[i][q] + alphas2[j][q] for q in range(n_vars))
                row.append(self.y[self.alpha_to_idx[summed]])
            rows.append(cp.hstack(row))
        M2 = cp.vstack(rows)
        constraints.append(M2 >> 0)

        # f-side bounds
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

        # SOC: v_k >= (a_k^2 + b_k^2)^2 = (y_{2e_k^a} + y_{2e_k^b})^2
        for k in range(K):
            ak2 = [0] * n_vars; ak2[k] = 2
            bk2 = [0] * n_vars; bk2[K + k] = 2
            sum_sq = self.y[self.alpha_to_idx[tuple(ak2)]] + self.y[self.alpha_to_idx[tuple(bk2)]]
            constraints.append(self.v[k] >= cp.square(sum_sq))

        # Fejer autocorrelation constraint
        weights = np.array([1.0 - k / (K + 1) for k in range(1, K + 1)])
        constraints.append(self.Omega >= 1.0 + 2.0 * (weights @ self.v))

        self.constraints = constraints
        self.problem = cp.Problem(cp.Minimize(self.Omega), constraints)

    def solve(self, solver: str = "MOSEK", verbose: bool = False, **kwargs) -> V20Result:
        val = self.problem.solve(solver=solver, verbose=verbose, **kwargs)
        return V20Result(
            status=self.problem.status,
            Omega=float(self.Omega.value) if self.Omega.value is not None else float("nan"),
            primal=float(val) if val is not None else float("nan"),
        )


if __name__ == "__main__":
    for N in (8, 16):
        for K in (3, 4):
            try:
                prob = AutocorrLowerBoundV20(N=N, K=K)
                out = prob.solve(solver="MOSEK", verbose=False)
                print(f"N={N:2d}, K={K:2d}: Omega={out.Omega:.6f}")
            except Exception as e:
                print(f"N={N}, K={K}: failed: {e}")
