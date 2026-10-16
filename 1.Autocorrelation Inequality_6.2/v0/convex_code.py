"""v0 baseline convex relaxation for the autocorrelation inequality (6.2).

Goal: prove a lower bound on the sharp constant

    C_{6.2} = inf_{f >= 0}  max_{|t|<=1/2} (f*f)(t) / (∫_{-1/4}^{1/4} f)^2.

Reduction.  Adding mass to f outside [-1/4, 1/4] does not change the
denominator (∫_{-1/4}^{1/4} f) while only increasing every value of f*f, so
the infimum is attained on functions supported in [-1/4, 1/4].  Homogeneity in
f lets us further normalize ∫f = 1.  Throughout this file we therefore work
with the normalized problem

    f >= 0,  supp(f) ⊂ [-1/4, 1/4],  ∫f = 1,
    F(t) := (f*f)(t),
    Ω    := max_{|t|<=1/2} F(t),

and we minimize Ω.  Any feasible value Ω of the relaxation is a valid lower
bound on C_{6.2}.

v0 SDP relaxation
-----------------

We discretize [-1/4, 1/4] into 2N intervals I_j of width L = 1/(4N) and let
p_j = ∫_{I_j} f.  Then p_j ≥ 0 and ∑p_j = 1.

f*f is bilinear in f, so we lift the rank-one outer product p p^T to a moment
matrix variable Q ∈ S^{2N×2N}.  The Schur-complement constraint

    [[1, p^T], [p, Q]]  ⪰  0

enforces Q ⪰ p p^T (i.e. every true f gives a feasible (p, Q = p p^T)).

For any band B_m = [-1/2 + (m-1)L, -1/2 + mL] of [-1/2, 1/2], with
m = 1, …, 4N, the integral of f*f over B_m is exactly bilinear in p:

    ∫_{B_m} (f*f)(t) dt
        = ∑_{j,k} (p_j/L)(p_k/L) · ∫_{B_m} (1_{I_j} * 1_{I_k})(t) dt
        = (1/L²) ∑_{j,k} p_j p_k · K_{m,j,k}

where K_{m,j,k} := ∫_{B_m}(1_{I_j} * 1_{I_k}).  For our equispaced grid the
convolution 1_{I_j}*1_{I_k} is a triangle of base 2L, height L, supported on
[-1/2 + (j+k-2)L, -1/2 + (j+k)L], split exactly in half by the boundary
between B_{j+k-1} and B_{j+k}, so

    K_{m,j,k} = L²/2  if  m ∈ {j+k-1, j+k},  else 0.

Substituting back, the constraint ∫_{B_m}(f*f) ≤ L · Ω becomes a linear
constraint in Q:

    (1/2) [ ∑_{j+k=m}  Q[j,k] + ∑_{j+k=m+1}  Q[j,k] ]  ≤  L · Ω.        (★)

The optimization problem is

    minimize     Ω
    subject to   p ≥ 0,  ∑ p_j = 1,
                 [[1, p^T], [p, Q]] ⪰ 0,
                 (★) for m = 1, …, 4N.

Validity (lower bound).  Every admissible f gives a feasible point with
Q = p p^T and Ω = ‖f*f‖_∞, since
    ∫_{B_m} f*f  ≤  |B_m| · ‖f*f‖_∞ = L · Ω.
Therefore the optimum of the relaxation is ≤ true infimum.  Any solver
optimum (or dual feasible certificate) is a valid lower bound on C_{6.2}.

Trivial value of v0.  Summing (★) over m and using the fact that each (j,k)
appears in exactly two band sums, we recover ∑_{j,k} Q[j,k] ≤ Ω.  Combined
with the PSD lift Q ⪰ p p^T and ∑p = 1, this gives Ω ≥ 1, the average
bound.  v0 therefore reproduces (and is meant to reproduce) the trivial
lower bound Ω* = 1.  Subsequent versions add new constraints to lift it.
"""

from __future__ import annotations

from dataclasses import dataclass

import cvxpy as cp
import numpy as np


@dataclass
class V0Result:
    status: str
    Omega: float
    primal: float


class AutocorrLowerBoundV0:
    def __init__(self, N: int = 8) -> None:
        if N < 1:
            raise ValueError("N must be a positive integer")
        self.N = N
        self.dim = 2 * N           # number of f-intervals on [-1/4, 1/4]
        self.L = 1.0 / (4 * N)     # interval width

        # decision variables
        self.Omega = cp.Variable(nonneg=True, name="Omega")
        self.p = cp.Variable(self.dim, nonneg=True, name="p")
        self.Q = cp.Variable((self.dim, self.dim), symmetric=True, name="Q")

        constraints: list[cp.constraints.Constraint] = []

        # PSD lift: [[1, p^T], [p, Q]] ⪰ 0
        block = cp.bmat([
            [np.array([[1.0]]), cp.reshape(self.p, (1, self.dim))],
            [cp.reshape(self.p, (self.dim, 1)), self.Q],
        ])
        constraints.append(block >> 0)

        # mass constraint
        constraints.append(cp.sum(self.p) == 1)

        # band constraints (★) for m = 1, ..., 4N
        n_bands = 4 * N
        # for efficiency, precompute lists of (j-1, k-1) indices for each j+k value
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

        self.constraints = constraints
        self.problem = cp.Problem(cp.Minimize(self.Omega), constraints)

    def solve(self, solver: str = "MOSEK", verbose: bool = False, **kwargs) -> V0Result:
        val = self.problem.solve(solver=solver, verbose=verbose, **kwargs)
        return V0Result(
            status=self.problem.status,
            Omega=float(self.Omega.value) if self.Omega.value is not None else float("nan"),
            primal=float(val) if val is not None else float("nan"),
        )


if __name__ == "__main__":
    prob = AutocorrLowerBoundV0(N=8)
    try:
        out = prob.solve(solver="MOSEK", verbose=True)
    except Exception as exc:
        print("MOSEK failed:", exc)
        out = prob.solve(solver="SCS", verbose=True)
    print()
    print(f"status      : {out.status}")
    print(f"primal value: {out.primal:.10f}")
    print(f"Omega       : {out.Omega:.10f}")
    print()
    print(f"v0 lower bound on C_6.2 : {out.Omega:.6f}")
