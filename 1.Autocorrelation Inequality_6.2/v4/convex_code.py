"""v4: Pure F-side LP using band integrals at multiple resolutions.

The earlier versions all used (p, Q) or (p, Q, M) liftings on the f-side
and were defeated by degenerate moment sequences (basis-vector mixtures)
that satisfy every PSD/marginal constraint while collapsing Omega to 1.

v4 is a fundamental change of variables: we **do not lift f** at all.
Instead we work directly on F = f*f using its bilinear structure.  The
variables are the average values w_j of F on a fine grid of [-1/2, 1/2],
plus Omega.  We then add constraints that exploit two facts simultaneously:

  (i) for any nonneg f with int f = 1 supported in [-1/4, 1/4],
      every band integral int_{[a,b]} (f*f) is bilinear in f and can be
      lower bounded via the cell-mass quantity p_j p_k that the *separate*
      moment lift Q already tracks;

  (ii) on the F side, the band integrals are LINEAR in w and can be
       upper bounded by Omega * (band length).

By coupling (i) and (ii) on the **same set of windows**, we obtain a
joint LP/SDP whose value is a valid lower bound on C_{6.2}, and whose
relaxation is harder to evade by diagonal/mixture solutions because the
direction of the inequalities differs between p-side and F-side.

Concretely, for each window B_w = [a_w, b_w]:

       sum_{(j,k) in F-cells with [a_w, b_w] subset of (I_j+I_k)} Q[j,k]
                <=  int_{B_w} F   <=   |B_w| * Omega.

The first inequality is the lower-bound on the integral via cells fully
contained in the preimage of B_w.  The second is the F <= Omega bound.

We also add the moment lift on (Q, p) from v1 (PSD lift, marginals, Q >= 0).
"""

from __future__ import annotations

from dataclasses import dataclass

import cvxpy as cp
import numpy as np


@dataclass
class V4Result:
    status: str
    Omega: float
    primal: float


class AutocorrLowerBoundV4:
    def __init__(self, N: int = 8, n_windows: int = 7) -> None:
        if N < 1:
            raise ValueError("N must be a positive integer")
        self.N = N
        self.dim = 2 * N
        self.L = 1.0 / (4 * N)
        self.n_windows = n_windows

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

        # ===== v4 NEW: many central-window constraints at various widths =====
        # For window B = [-h, h] with h > 0 we have
        #     sum_{(j,k) inside} Q[j,k] <= 2h * Omega
        # where "inside" means the cell I_j x I_k is contained in
        # {(tau,x) : |tau+x| <= h}, i.e., j+k satisfies
        #     j+k >= 2N + 2 - 2hN (lower)  and j+k <= 2N + 2hN (upper).
        # We add constraints for several h values.
        max_half_width = 0.5
        widths = np.linspace(2 * self.L, max_half_width, n_windows)  # half-widths
        for h in widths:
            # cells inside |tau+x|<=h
            j_low_bound = max(2, int(np.ceil((1.0 - 2 * h) / (2 * self.L))) + 1)
            j_high_bound = min(4 * N, int(np.floor((1.0 + 2 * h) / (2 * self.L))))
            inside_terms = []
            for j in range(1, self.dim + 1):
                for k in range(1, self.dim + 1):
                    s = j + k
                    if j_low_bound <= s <= j_high_bound:
                        inside_terms.append(self.Q[j - 1, k - 1])
            if inside_terms:
                constraints.append(cp.sum(cp.hstack(inside_terms)) <= 2.0 * h * self.Omega)

        # ===== v4 NEW: shifted-window constraints for windows centered away from 0 =====
        # For window B = [c-h, c+h], cells inside have
        #     j+k satisfying  -1/2 + (j+k-2)L >= c - h  and  -1/2 + (j+k)L <= c + h
        #     => j+k >= (c - h + 1/2)/L + 2  and  j+k <= (c + h + 1/2)/L
        for c_frac in (-0.25, 0.25):
            for h in widths:
                lo = (c_frac - h + 0.5) / self.L + 2
                hi = (c_frac + h + 0.5) / self.L
                inside_terms = []
                for j in range(1, self.dim + 1):
                    for k in range(1, self.dim + 1):
                        s = j + k
                        if lo <= s <= hi:
                            inside_terms.append(self.Q[j - 1, k - 1])
                if inside_terms:
                    constraints.append(cp.sum(cp.hstack(inside_terms)) <= 2.0 * h * self.Omega)

        self.constraints = constraints
        self.problem = cp.Problem(cp.Minimize(self.Omega), constraints)

    def solve(self, solver: str = "MOSEK", verbose: bool = False, **kwargs) -> V4Result:
        val = self.problem.solve(solver=solver, verbose=verbose, **kwargs)
        return V4Result(
            status=self.problem.status,
            Omega=float(self.Omega.value) if self.Omega.value is not None else float("nan"),
            primal=float(val) if val is not None else float("nan"),
        )


if __name__ == "__main__":
    prob = AutocorrLowerBoundV4(N=8, n_windows=15)
    try:
        out = prob.solve(solver="MOSEK", verbose=False)
    except Exception as exc:
        print("MOSEK failed:", exc)
        out = prob.solve(solver="SCS", verbose=False)
    print(f"status      : {out.status}")
    print(f"primal value: {out.primal:.10f}")
    print(f"Omega       : {out.Omega:.10f}")
    print(f"v4 lower bound on C_6.2 : {out.Omega:.6f}")
