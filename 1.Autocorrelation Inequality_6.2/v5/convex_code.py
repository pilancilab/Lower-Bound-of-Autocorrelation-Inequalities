"""v5 convex relaxation: even-f restricted with anti-diagonal-equals-diagonal constraint.

Key insight that breaks the diagonal-Q degeneracy
--------------------------------------------------

The previous relaxations all admitted the pathological feasible point
   p_j = 1/(2N),   Q[j,k] = (1/(2N)) * delta_{jk}    (j,k = 1..2N)
with Omega = 1.  This Q has every diagonal entry equal to 1/(2N) and every
off-diagonal entry equal to 0.

For an admissible f that is **even** (i.e. f(-x) = f(x)) and piecewise
constant on the symmetric grid I_j, j=1..2N, we have p_j = p_{2N+1-j}
(masses are reflection-symmetric).  Then for the rank-one Q = p p^T:

   Q[j, 2N+1-j]  =  p_j p_{2N+1-j}  =  p_j * p_j  =  p_j^2  =  Q[j, j].

So the **anti-diagonal of Q (with respect to the cell-reflection
(j, 2N+1-j)) coincides with the diagonal of Q** for any rank-one Q from
an even f.

Now check the diagonal Q above:
   Q[j, 2N+1-j]  =  (1/(2N)) * delta_{j, 2N+1-j}  =  0     (since 2j != 2N+1)
   Q[j, j]       =  1/(2N).
The two are not equal, so the constraint Q[j, 2N+1-j] = Q[j, j] **kills the
degenerate diagonal point** while remaining valid for every rank-one Q
coming from an admissible even f.

What v5 proves
--------------

Adding this constraint, plus
  (i) the even-f symmetry constraint p_j = p_{2N+1-j},
  (ii) the relation Q[j, k] = Q[2N+1-j, 2N+1-k] (Q is reflection-symmetric),
yields a strictly tighter relaxation whose optimum lower-bounds
**the even-restricted constant**

   C^{even}_{6.2}  :=  inf {  ||f * f||_infty  :  f admissible AND f even  }.

Caveat: C^{even}_{6.2} >= C_{6.2} (since restricting f's gives a larger inf).
The bound from v5 is therefore a lower bound on C^{even}_{6.2}, NOT
directly on C_{6.2}.  However the relaxation as a "sanity check" verifies
that the analytic bound C^{even}_{6.2} >= 2 (proved via Cauchy-Schwarz
on f) can be (re)proved by the convex program.
"""

from __future__ import annotations

from dataclasses import dataclass

import cvxpy as cp
import numpy as np


@dataclass
class V5Result:
    status: str
    Omega: float
    primal: float


class AutocorrLowerBoundV5:
    def __init__(self, N: int = 8) -> None:
        if N < 1:
            raise ValueError("N must be a positive integer")
        self.N = N
        self.dim = 2 * N
        self.L = 1.0 / (4 * N)

        self.Omega = cp.Variable(nonneg=True, name="Omega")
        self.p = cp.Variable(self.dim, nonneg=True, name="p")
        self.Q = cp.Variable((self.dim, self.dim), symmetric=True, name="Q")

        constraints: list = []

        # PSD lift
        block = cp.bmat([
            [np.array([[1.0]]), cp.reshape(self.p, (1, self.dim))],
            [cp.reshape(self.p, (self.dim, 1)), self.Q],
        ])
        constraints.append(block >> 0)

        constraints.append(cp.sum(self.p) == 1)
        constraints.append(self.Q >= 0)
        constraints.append(cp.sum(self.Q, axis=1) == self.p)

        # ===== v5 NEW: even-f symmetry =====
        # p_j = p_{2N+1-j}    (j-th vs (2N+1-j)-th interval in 1-indexed; in 0-indexed
        # this is p[i] = p[dim - 1 - i] for i = 0..dim-1)
        for i in range(self.dim // 2):
            constraints.append(self.p[i] == self.p[self.dim - 1 - i])

        # Q[j, k] = Q[2N+1-j, 2N+1-k] (in 0-indexed: Q[i, l] = Q[dim-1-i, dim-1-l])
        for i in range(self.dim):
            for l in range(self.dim):
                if i + l < self.dim - 1:
                    # add only one direction; symmetry handles the other
                    i2 = self.dim - 1 - i
                    l2 = self.dim - 1 - l
                    if (i2, l2) != (i, l):
                        constraints.append(self.Q[i, l] == self.Q[i2, l2])

        # ===== v5 NEW: anti-diagonal-equals-diagonal =====
        # Q[j, 2N+1-j] = Q[j, j]    (in 0-indexed: Q[i, dim-1-i] = Q[i, i])
        for i in range(self.dim):
            constraints.append(self.Q[i, self.dim - 1 - i] == self.Q[i, i])

        # base v0 band constraints
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

        # central window from v1
        center_terms = []
        for j in range(1, self.dim + 1):
            for k in range(1, self.dim + 1):
                if N + 2 <= j + k <= 3 * N:
                    center_terms.append(self.Q[j - 1, k - 1])
        if center_terms:
            constraints.append(cp.sum(cp.hstack(center_terms)) <= 0.5 * self.Omega)

        self.constraints = constraints
        self.problem = cp.Problem(cp.Minimize(self.Omega), constraints)

    def solve(self, solver: str = "MOSEK", verbose: bool = False, **kwargs) -> V5Result:
        val = self.problem.solve(solver=solver, verbose=verbose, **kwargs)
        return V5Result(
            status=self.problem.status,
            Omega=float(self.Omega.value) if self.Omega.value is not None else float("nan"),
            primal=float(val) if val is not None else float("nan"),
        )


if __name__ == "__main__":
    prob = AutocorrLowerBoundV5(N=6)
    try:
        out = prob.solve(solver="MOSEK", verbose=False)
    except Exception as exc:
        print("MOSEK failed:", exc)
        out = prob.solve(solver="SCS", verbose=False)
    print(f"status      : {out.status}")
    print(f"primal value: {out.primal:.10f}")
    print(f"Omega       : {out.Omega:.10f}")
    print(f"v5 lower bound on C^{{even}}_6.2 : {out.Omega:.6f}")
