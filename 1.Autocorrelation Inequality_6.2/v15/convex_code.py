"""v15: v14 autocorrelation Fejer combined with v7/v11 windowed (W) constraints.

Builds on v14's rigorous autocorrelation Fejer test (valid for all
admissible f, not just even) and adds:

1. The (p, Q) SDP lift from v0/v7 (p_j = mass on j-th f-cell, Q[j,k]
   PSD-lifted as a moment matrix).
2. The rigorous v7 (W) windowed integral inequalities at multiple
   widths and offsets.

The combination is rigorous because all v14 and v7/v11 constraints
are individually valid for any admissible f, and adding them only
restricts the relaxation feasible set without excluding any
admissible point.

The expected effect: v14 alone gives Omega ~ 1.06 at small (N, K).
The (W) constraints constrain Q (which encodes p_j p_k), and via the
linear bounds on a_k, b_k from p, this propagates to tighter bounds
on a_k^2 + b_k^2 = u_k.  The Fejer constraint then gives a larger
Omega.
"""

from __future__ import annotations

from dataclasses import dataclass

import cvxpy as cp
import numpy as np


@dataclass
class V15Result:
    status: str
    Omega: float
    primal: float


class AutocorrLowerBoundV15:
    def __init__(self, N: int = 8, K: int = 8) -> None:
        if N < 1 or K < 1:
            raise ValueError("N, K must be positive")
        self.N = N
        self.K = K
        self.dim = 2 * N
        self.L = 1.0 / (4 * N)

        self.Omega = cp.Variable(nonneg=True, name="Omega")
        self.p = cp.Variable(self.dim, nonneg=True, name="p")
        self.Q = cp.Variable((self.dim, self.dim), symmetric=True, name="Q")
        self.a = cp.Variable(K + 1, name="a")
        self.b = cp.Variable(K + 1, name="b")
        self.M = cp.Variable((2 * K + 1, 2 * K + 1), symmetric=True, name="M")
        self.v = cp.Variable(K, nonneg=True, name="v")

        constraints: list = []

        # ---- (p, Q) PSD lift + standard ----
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

        # ---- M PSD lift on (1, a, b) ----
        constraints.append(self.M >> 0)
        constraints.append(self.M[0, 0] == 1)
        for k in range(1, K + 1):
            constraints.append(self.M[0, k] == self.a[k])
            constraints.append(self.M[k, 0] == self.a[k])
            constraints.append(self.M[0, K + k] == self.b[k])
            constraints.append(self.M[K + k, 0] == self.b[k])

        # ---- f-side bounds on a_k, b_k from p ----
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

        # ---- v14 autocorrelation Fejer constraint ----
        for k in range(1, K + 1):
            sum_diag = self.M[k, k] + self.M[K + k, K + k]
            constraints.append(self.v[k - 1] >= cp.square(sum_diag))
        weights = np.array([1.0 - k / (K + 1) for k in range(1, K + 1)])
        constraints.append(self.Omega >= 1.0 + 2.0 * (weights @ self.v))

        # ---- v7/v11 (W) windowed constraints (centered + shifted) ----
        widths = list(np.linspace(2.0 * self.L, 0.5, max(8, self.dim)))
        offsets = [0.0] + [c * self.L for c in range(1, self.dim // 2 + 1)] + [-c * self.L for c in range(1, self.dim // 2 + 1)]
        for c in offsets:
            for h in widths:
                a_w = c - h
                b_w = c + h
                if a_w < -0.5 or b_w > 0.5:
                    continue
                lo = int(np.ceil((a_w + 0.5) / self.L)) + 2
                hi = int(np.floor((b_w + 0.5) / self.L))
                inside = []
                for j in range(1, self.dim + 1):
                    for k in range(1, self.dim + 1):
                        if lo <= j + k <= hi:
                            inside.append(self.Q[j - 1, k - 1])
                if inside:
                    constraints.append(cp.sum(cp.hstack(inside)) <= 2.0 * h * self.Omega)

        self.constraints = constraints
        self.problem = cp.Problem(cp.Minimize(self.Omega), constraints)

    def solve(self, solver: str = "MOSEK", verbose: bool = False, **kwargs) -> V15Result:
        val = self.problem.solve(solver=solver, verbose=verbose, **kwargs)
        return V15Result(
            status=self.problem.status,
            Omega=float(self.Omega.value) if self.Omega.value is not None else float("nan"),
            primal=float(val) if val is not None else float("nan"),
        )


if __name__ == "__main__":
    for N in (8, 12):
        for K in (8, 16):
            try:
                prob = AutocorrLowerBoundV15(N=N, K=K)
                out = prob.solve(solver="MOSEK", verbose=False)
                print(f"N={N:2d}, K={K:2d}: Omega={out.Omega:.6f}")
            except Exception as e:
                print(f"N={N}, K={K}: failed: {e}")
