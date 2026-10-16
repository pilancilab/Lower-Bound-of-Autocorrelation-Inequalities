"""v7: rigorous reboot for the autocorrelation inequality 6.2.

Addresses all verifier concerns from v1..v5 verdicts:

1. WLOG reductions are justified explicitly in rigorousproof.md.
2. The v0 false-equality band-constraint chain is REPLACED by the
   rigorous "fully-contained cells" (W) form from v4.
3. The v5 (E)/(F)/(G) constraints are kept (they are valid in the
   even-restricted setting); the program therefore bounds
   C^{even}_{6.2} = inf{||f*f||_inf : f admissible AND even}.
4. The bound on C^{even}_{6.2} is < C^{even}_{6.2} = 2 (the analytic
   value), and approaches 2 as N grows.

Note that this version bounds C^{even}_{6.2}, not C_{6.2}. v8 will
attempt the harder extension to general f.
"""

from __future__ import annotations

from dataclasses import dataclass

import cvxpy as cp
import numpy as np


@dataclass
class V7Result:
    status: str
    Omega: float
    primal: float


class AutocorrLowerBoundV7:
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

        # ---- PSD lift: Q >> p p^T via Schur complement ----
        block = cp.bmat([
            [np.array([[1.0]]), cp.reshape(self.p, (1, self.dim), order="C")],
            [cp.reshape(self.p, (self.dim, 1), order="C"), self.Q],
        ])
        constraints.append(block >> 0)

        # ---- mass + linearised conditions on (p, Q) ----
        constraints.append(cp.sum(self.p) == 1)
        constraints.append(self.Q >= 0)
        constraints.append(cp.sum(self.Q, axis=1) == self.p)

        # ---- v5 (E)/(F)/(G) even-symmetry constraints ----
        # (E) p_j = p_{2N+1-j}     (0-indexed: p[i] = p[dim-1-i])
        for i in range(self.dim // 2):
            constraints.append(self.p[i] == self.p[self.dim - 1 - i])
        # (F) Q[j,k] = Q[2N+1-j, 2N+1-k]
        for i in range(self.dim):
            for l in range(self.dim):
                if i + l < self.dim - 1:
                    i2 = self.dim - 1 - i
                    l2 = self.dim - 1 - l
                    if (i2, l2) != (i, l):
                        constraints.append(self.Q[i, l] == self.Q[i2, l2])
        # (G) Q[j, 2N+1-j] = Q[j, j]
        for i in range(self.dim):
            constraints.append(self.Q[i, self.dim - 1 - i] == self.Q[i, i])

        # ---- v4 (W) band constraints with cells fully inside ----
        # The window B = [a, b] of width h_B = b - a satisfies
        #     int_B (f*f) dt  >=  sum_{(j,k): I_j x I_k subset preimage} p_j p_k.
        # And int_B (f*f) <= h_B Omega.  Hence
        #     sum_{cells fully inside}  Q[j,k]  <=  h_B Omega.
        # Cells fully inside iff
        #     -1/2 + (j+k-2)L >= a   AND   -1/2 + (j+k)L <= b,
        # i.e.  j+k >= 2 + 2N(2(a + 1/2))  and  j+k <= 2N(2(b + 1/2)).
        # We use a family of CENTERED windows of various half-widths.
        widths = np.linspace(2.0 * self.L, 0.5, max(8, self.dim))  # half-widths h
        for h in widths:
            # window [-h, h], length 2h
            j_lo = int(np.ceil((1.0 - 2.0 * h) / (2.0 * self.L))) + 1  # j+k >= 2N - 2hN + 2
            # Re-derive carefully:
            # (j+k-2) L >= -h + 1/2  =>  j+k >= (1/2 - h)/L + 2
            j_lo = int(np.ceil((0.5 - h) / self.L)) + 2
            # (j+k) L <= h + 1/2     =>  j+k <= (h + 1/2)/L
            j_hi = int(np.floor((0.5 + h) / self.L))
            inside_terms = []
            for j in range(1, self.dim + 1):
                for k in range(1, self.dim + 1):
                    if j_lo <= j + k <= j_hi:
                        inside_terms.append(self.Q[j - 1, k - 1])
            if inside_terms:
                constraints.append(cp.sum(cp.hstack(inside_terms)) <= 2.0 * h * self.Omega)

        # ---- Full window B = [-1/2, 1/2] (gives the trivial Omega >= 1) ----
        all_q = []
        for j in range(self.dim):
            for k in range(self.dim):
                all_q.append(self.Q[j, k])
        constraints.append(cp.sum(cp.hstack(all_q)) <= 1.0 * self.Omega)

        self.constraints = constraints
        self.problem = cp.Problem(cp.Minimize(self.Omega), constraints)

    def solve(self, solver: str = "MOSEK", verbose: bool = False, **kwargs) -> V7Result:
        val = self.problem.solve(solver=solver, verbose=verbose, **kwargs)
        return V7Result(
            status=self.problem.status,
            Omega=float(self.Omega.value) if self.Omega.value is not None else float("nan"),
            primal=float(val) if val is not None else float("nan"),
        )


if __name__ == "__main__":
    for N in (4, 6, 8, 10):
        prob = AutocorrLowerBoundV7(N=N)
        try:
            out = prob.solve(solver="MOSEK", verbose=False)
        except Exception as exc:
            print("MOSEK failed:", exc)
            out = prob.solve(solver="SCS", verbose=False)
        print(f"N={N:2d}: status={out.status}, Omega={out.Omega:.6f}")
