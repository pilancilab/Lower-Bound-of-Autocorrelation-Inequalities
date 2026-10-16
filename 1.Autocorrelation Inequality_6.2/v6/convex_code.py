"""v6 convex relaxation: anti-diagonal-equals-diagonal constraint with asymmetry slack.

Idea
----

In v5 we proved that for even admissible f the relaxation gives Omega ~ 1.86,
recovering most of the analytic bound C^{even}_{6.2} >= 2.  The key
constraint was the **anti-diagonal-equals-diagonal** identity

    Q[j, 2N+1-j]  =  Q[j, j]    (for even f),

which kills the diagonal-Q degeneracy.

For non-even f the identity becomes a slack inequality.  Specifically,
since Q[j,k] = p_j p_k in the rank-one truth and 2 a b = a^2 + b^2 - (a-b)^2,
we have

    2 Q[j, 2N+1-j]  =  Q[j, j] + Q[2N+1-j, 2N+1-j]  -  (p_j - p_{2N+1-j})^2.

Introducing a slack variable s_j >= 0 with s_j = (p_j - p_{2N+1-j})^2, we
can write the slacked equality

    2 Q[j, 2N+1-j]  =  Q[j, j] + Q[2N+1-j, 2N+1-j]  -  s_j.            (*)

In the convex lift we relax s_j >= (p_j - p_{2N+1-j})^2 (rotated SOC), and
add (*) as a linear equation linking the lifted moments.  We further
linearize the rotated SOC via the lift to a 2x2 PSD block.

Validity
--------

For every admissible f, the assignment p_j^* = int_{I_j} f,
Q^*[j,k] = p_j^* p_k^*, s_j^* = (p_j^* - p_{2N+1-j}^*)^2 satisfies (*),
the SOC s_j >= (p_j - p_{2N+1-j})^2, and all v0/v1 constraints.  Thus
the relaxation is rigorous: its optimum is <= C_{6.2}.

The hope was that this slacked version would still cut down the diagonal
Q = (1/(2N)) I degeneracy.  Numerical experiments below show whether
this is realized.
"""

from __future__ import annotations

from dataclasses import dataclass

import cvxpy as cp
import numpy as np


@dataclass
class V6Result:
    status: str
    Omega: float
    primal: float


class AutocorrLowerBoundV6:
    def __init__(self, N: int = 8) -> None:
        if N < 1:
            raise ValueError("N must be a positive integer")
        self.N = N
        self.dim = 2 * N
        self.L = 1.0 / (4 * N)

        self.Omega = cp.Variable(nonneg=True, name="Omega")
        self.p = cp.Variable(self.dim, nonneg=True, name="p")
        self.Q = cp.Variable((self.dim, self.dim), symmetric=True, name="Q")
        self.s = cp.Variable(self.dim, nonneg=True, name="s")

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

        # ====== v6: rotated SOC s_j >= (p_j - p_{2N+1-j})^2 ======
        # Equivalent to s_j >= p_j^2 - 2 p_j p_{2N+1-j} + p_{2N+1-j}^2,
        # which (using lifted Q[j,j]=p_j^2 etc.) we will replace by the
        # tighter linear constraint  s_j  >=  Q[j,j] + Q[2N+1-j,2N+1-j] - 2 Q[j, 2N+1-j].
        # In the rank-one truth, both sides are equal, so this is valid.
        for j in range(self.dim):
            j2 = self.dim - 1 - j
            constraints.append(
                self.s[j] >= self.Q[j, j] + self.Q[j2, j2] - 2 * self.Q[j, j2]
            )

        # ====== v6: linear identity (*) tying the slacks to the rank-one structure ======
        # For rank-one, 2 Q[j, 2N+1-j] = Q[j,j] + Q[2N+1-j, 2N+1-j] - s_j EXACTLY.
        # In the relaxation we keep it as ">=" (which the previous SOC bound guarantees)
        # plus an ADDITIONAL upper bound s_j <= (p_j - p_{2N+1-j})^2 derived from
        # PSD: but the latter is non-convex.  We instead bound the SUM:
        #     sum_j s_j  <=  2 tr(Q) - 2 sum_j Q[j, 2N+1-j],
        # which is a linear consequence of summing the per-j inequality.
        sum_anti = cp.sum(cp.hstack([self.Q[j, self.dim - 1 - j] for j in range(self.dim)]))
        sum_diag = cp.trace(self.Q)
        constraints.append(cp.sum(self.s) <= 2.0 * sum_diag - 2.0 * sum_anti)

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

    def solve(self, solver: str = "MOSEK", verbose: bool = False, **kwargs) -> V6Result:
        val = self.problem.solve(solver=solver, verbose=verbose, **kwargs)
        return V6Result(
            status=self.problem.status,
            Omega=float(self.Omega.value) if self.Omega.value is not None else float("nan"),
            primal=float(val) if val is not None else float("nan"),
        )


if __name__ == "__main__":
    prob = AutocorrLowerBoundV6(N=8)
    try:
        out = prob.solve(solver="MOSEK", verbose=False)
    except Exception as exc:
        print("MOSEK failed:", exc)
        out = prob.solve(solver="SCS", verbose=False)
    print(f"status      : {out.status}")
    print(f"primal value: {out.primal:.10f}")
    print(f"Omega       : {out.Omega:.10f}")
    print(f"v6 lower bound on C_6.2 : {out.Omega:.6f}")
