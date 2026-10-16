"""v1 convex relaxation for the autocorrelation inequality (6.2).

Improvement over v0
-------------------

(A) Element-wise nonnegativity:  Q[j,k] >= 0 for all j, k.
    Justification: p_j p_k >= 0 in the rank-one truth.

(B) Row-sum (marginal) constraint:  sum_k Q[j,k] = p_j  for all j.
    Justification: with sum_k p_k = 1, sum_k p_j p_k = p_j.

(C) Central-window integral constraint (KEY):

        sum_{(j,k):  I_j x I_k subset of {(tau,x): |tau+x| <= 1/4}}  Q[j,k]
            <=  Omega / 2.

    Justification: integrating (f*f) over the central window [-1/4, 1/4]
    gives a quantity bounded above by Omega · (length of window) = Omega/2,
    and bounded below by the cell-mass of every (j,k) whose product cell is
    entirely contained in the preimage of the window.  See rigorousproof.md.

The rest of the program is identical to v0.
"""

from __future__ import annotations

from dataclasses import dataclass

import cvxpy as cp
import numpy as np


@dataclass
class V1Result:
    status: str
    Omega: float
    primal: float


class AutocorrLowerBoundV1:
    def __init__(self, N: int = 8) -> None:
        if N < 1:
            raise ValueError("N must be a positive integer")
        self.N = N
        self.dim = 2 * N
        self.L = 1.0 / (4 * N)

        self.Omega = cp.Variable(nonneg=True, name="Omega")
        self.p = cp.Variable(self.dim, nonneg=True, name="p")
        self.Q = cp.Variable((self.dim, self.dim), symmetric=True, name="Q")

        constraints: list[cp.constraints.Constraint] = []

        # PSD lift
        block = cp.bmat([
            [np.array([[1.0]]), cp.reshape(self.p, (1, self.dim))],
            [cp.reshape(self.p, (self.dim, 1)), self.Q],
        ])
        constraints.append(block >> 0)

        # mass
        constraints.append(cp.sum(self.p) == 1)

        # (A) nonneg, (B) marginals
        constraints.append(self.Q >= 0)
        constraints.append(cp.sum(self.Q, axis=1) == self.p)

        # band constraints
        n_bands = 4 * N
        index_lists: dict[int, list[tuple[int, int]]] = {}
        for j in range(1, self.dim + 1):
            for k in range(1, self.dim + 1):
                index_lists.setdefault(j + k, []).append((j - 1, k - 1))

        for m in range(1, n_bands + 1):
            terms = []
            for s in (m, m + 1):
                if s in index_lists:
                    for (jj, kk) in index_lists[s]:
                        terms.append(self.Q[jj, kk])
            if terms:
                constraints.append(0.5 * cp.sum(cp.hstack(terms)) <= self.L * self.Omega)

        # (C) central-window integral constraint:
        #   sum over (j,k) with cell I_j x I_k contained in {|tau+x| <= 1/4}
        #   of Q[j,k]  <=  Omega / 2.
        # The cell I_j x I_k is contained in {|tau+x| <= 1/4} iff
        #     min(tau+x) = -1/2 + (j+k-2)L >= -1/4   and
        #     max(tau+x) = -1/2 + (j+k)L   <=  1/4,
        # i.e.  N+2 <= j+k <= 3N.
        center_terms = []
        for j in range(1, self.dim + 1):
            for k in range(1, self.dim + 1):
                if N + 2 <= j + k <= 3 * N:
                    center_terms.append(self.Q[j - 1, k - 1])
        if center_terms:
            constraints.append(cp.sum(cp.hstack(center_terms)) <= 0.5 * self.Omega)

        self.constraints = constraints
        self.problem = cp.Problem(cp.Minimize(self.Omega), constraints)

    def solve(self, solver: str = "MOSEK", verbose: bool = False, **kwargs) -> V1Result:
        val = self.problem.solve(solver=solver, verbose=verbose, **kwargs)
        return V1Result(
            status=self.problem.status,
            Omega=float(self.Omega.value) if self.Omega.value is not None else float("nan"),
            primal=float(val) if val is not None else float("nan"),
        )


if __name__ == "__main__":
    prob = AutocorrLowerBoundV1(N=8)
    try:
        out = prob.solve(solver="MOSEK", verbose=False)
    except Exception as exc:
        print("MOSEK failed:", exc)
        out = prob.solve(solver="SCS", verbose=False)
    print(f"status      : {out.status}")
    print(f"primal value: {out.primal:.10f}")
    print(f"Omega       : {out.Omega:.10f}")
    print(f"v1 lower bound on C_6.2 : {out.Omega:.6f}")
