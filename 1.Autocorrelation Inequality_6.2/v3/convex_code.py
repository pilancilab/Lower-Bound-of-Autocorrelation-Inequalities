"""v3: Lasserre level-2 SOS relaxation for the autocorrelation inequality 6.2.

Why level 2 and not level 1?
----------------------------

v0/v1/v2 are level-1 (Shor / DNN) relaxations.  Their feasible set is
strictly larger than the convex hull of rank-one outer products of
nonneg vectors, which lets the relaxation choose pathological matrices
like Q = (1/(2N)) * I that satisfy every linear and PSD constraint while
having Omega = 1.  These pathological matrices do not correspond to any
true rank-one product p p^T, but level-1 cannot tell.

The Lasserre level-2 hierarchy adds:

  - All moments y_alpha for |alpha| <= 4 (in particular, third- and
    fourth-order moments of p).
  - The level-2 moment matrix M_2 of size {monomials of degree <= 2},
    constrained PSD.
  - Localizing matrices L_j for the constraints p_j >= 0, of size
    {monomials of degree <= 1}, constrained PSD.

These collectively force the joint moments of p to be consistent with
SOME nonneg vector p, which strictly tightens the relaxation.

Computational note: level 2 grows fast.  For n = 2N variables,
  - moments of degree <= 4: C(n+4, 4)
  - moment matrix size:    C(n+2, 2)
  - localizing matrix size: n+1, with n constraints
We pick a small N (default N = 4, n = 8) so that the SDP stays tractable.

Notation
--------

Let n = 2N.  We index multi-indices alpha = (a_1, ..., a_n) with |alpha|
the sum of components.  Moments y_alpha = E[ p^alpha ] in the truth, but
in the relaxation y_alpha is a free variable with PSD/linear constraints.

Conventions:
  y_0  := y_{(0,0,...,0)} = 1                         (mass)
  p_j  := y_{e_j}                                       (degree 1)
  Q[j,k] := y_{e_j + e_k}                               (degree 2)

The discretization of f on [-1/4, 1/4] uses 2N intervals of width
L = 1/(4N).  The same band-integral derivation as v0 gives the
autocorrelation constraints

   0.5 * sum_{(j,k): j+k in {m, m+1}}  Q[j,k]  <=  L * Omega,
   for m = 1, ..., 4N.

The validity proof is identical to v0; the level-2 lift only sharpens
which (Q, p) tuples are admissible.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations_with_replacement

import cvxpy as cp
import numpy as np


@dataclass
class V3Result:
    status: str
    Omega: float
    primal: float
    n_moments: int


def enumerate_alphas(n: int, max_deg: int) -> list[tuple[int, ...]]:
    out = []
    for d in range(max_deg + 1):
        for combo in combinations_with_replacement(range(n), d):
            alpha = [0] * n
            for i in combo:
                alpha[i] += 1
            out.append(tuple(alpha))
    return out


def add_alphas(a: tuple[int, ...], b: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(x + y for x, y in zip(a, b))


class AutocorrLowerBoundV3:
    def __init__(self, N: int = 4) -> None:
        if N < 1:
            raise ValueError("N must be a positive integer")
        self.N = N
        self.n = 2 * N
        self.L = 1.0 / (4 * N)

        # ---- enumerate multi-indices of degree <= 4 ----
        self.alphas4 = enumerate_alphas(self.n, 4)
        self.alpha_to_idx = {a: i for i, a in enumerate_alphas_with_idx(self.alphas4)}

        # Build a single CVXPY Variable indexed by alpha (one entry per moment).
        n_moments = len(self.alphas4)
        self.y = cp.Variable(n_moments, name="y")
        self.n_moments = n_moments

        # Convenience: index of the empty multi-index (degree 0)
        self.zero_alpha = tuple([0] * self.n)
        self.zero_idx = self.alpha_to_idx[self.zero_alpha]

        # Index of e_j (degree 1)
        self.e = []
        for j in range(self.n):
            ej = [0] * self.n
            ej[j] = 1
            self.e.append(tuple(ej))

        # ---- variables for Omega ----
        self.Omega = cp.Variable(nonneg=True, name="Omega")

        constraints: list = []

        # ---- normalization ----
        constraints.append(self.y[self.zero_idx] == 1)

        # ---- mass: sum_j p_j = 1 ----
        constraints.append(
            cp.sum(cp.hstack([self.y[self.alpha_to_idx[ej]] for ej in self.e])) == 1
        )

        # ---- y[alpha] are free except y_0 = 1 (sym constraint not needed: y is a vector) ----

        # ---- Build moment matrix M_2 of size |alphas with deg <= 2| ----
        alphas2 = [a for a in self.alphas4 if sum(a) <= 2]
        # M_2[i, j] = y_{alpha_i + alpha_j}
        size2 = len(alphas2)
        rows = []
        for i in range(size2):
            row = []
            for j in range(size2):
                s = add_alphas(alphas2[i], alphas2[j])
                row.append(self.y[self.alpha_to_idx[s]])
            rows.append(cp.hstack(row))
        M_2 = cp.vstack(rows)
        constraints.append(M_2 == M_2.T)  # symmetry (will be auto)
        constraints.append(M_2 >> 0)

        # ---- Localizing matrices for p_j >= 0 ----
        alphas1 = [a for a in self.alphas4 if sum(a) <= 1]
        size1 = len(alphas1)
        for jdx in range(self.n):
            ej = self.e[jdx]
            rows_j = []
            for i in range(size1):
                row = []
                for k in range(size1):
                    s = add_alphas(add_alphas(alphas1[i], alphas1[k]), ej)
                    row.append(self.y[self.alpha_to_idx[s]])
                rows_j.append(cp.hstack(row))
            L_j = cp.vstack(rows_j)
            constraints.append(L_j >> 0)

        # ---- Localizing matrices for sum_j p_j = 1 (equality) ----
        # For each alpha of degree <= 2: sum_j y_{alpha + e_j} = y_alpha
        for a in alphas2:
            constraints.append(
                cp.sum(cp.hstack([self.y[self.alpha_to_idx[add_alphas(a, ej)]] for ej in self.e])) == self.y[self.alpha_to_idx[a]]
            )

        # ---- Band autocorrelation constraints ----
        # Q[j,k] = y_{e_j + e_k}
        def Q(j: int, k: int):
            return self.y[self.alpha_to_idx[add_alphas(self.e[j], self.e[k])]]

        n_bands = 4 * N
        index_lists: dict[int, list[tuple[int, int]]] = {}
        for j in range(1, self.n + 1):
            for k in range(1, self.n + 1):
                index_lists.setdefault(j + k, []).append((j - 1, k - 1))
        for m in range(1, n_bands + 1):
            terms = []
            for s in (m, m + 1):
                if s in index_lists:
                    for (jj, kk) in index_lists[s]:
                        terms.append(Q(jj, kk))
            if terms:
                constraints.append(0.5 * cp.sum(cp.hstack(terms)) <= self.L * self.Omega)

        self.constraints = constraints
        self.problem = cp.Problem(cp.Minimize(self.Omega), constraints)

    def solve(self, solver: str = "MOSEK", verbose: bool = False, **kwargs) -> V3Result:
        val = self.problem.solve(solver=solver, verbose=verbose, **kwargs)
        return V3Result(
            status=self.problem.status,
            Omega=float(self.Omega.value) if self.Omega.value is not None else float("nan"),
            primal=float(val) if val is not None else float("nan"),
            n_moments=self.n_moments,
        )


def enumerate_alphas_with_idx(alphas):
    return [(i, a) for i, a in enumerate(alphas)]
# Note: helper above renamed inline; provided for clarity in alpha_to_idx setup.


if __name__ == "__main__":
    prob = AutocorrLowerBoundV3(N=3)
    print(f"# moments: {prob.n_moments}")
    try:
        out = prob.solve(solver="MOSEK", verbose=False)
    except Exception as exc:
        print("MOSEK failed:", exc)
        out = prob.solve(solver="SCS", verbose=False)
    print(f"status      : {out.status}")
    print(f"primal value: {out.primal:.10f}")
    print(f"Omega       : {out.Omega:.10f}")
    print(f"v3 lower bound on C_6.2 : {out.Omega:.6f}")
