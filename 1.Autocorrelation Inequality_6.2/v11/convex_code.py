"""v11: v9 Fejer-kernel test extended with shifted windows and v7 (W).

Per the v7 verdict's Comment 3 / Feedback Item 1, v9's bound on
C^{even}_{6.2} can be sharpened by adding *more windows* (shifted and
asymmetric) on top of the centered family.  v11 combines:

  - v7's rigorous (W) family of central windows (cells fully inside)
    plus shifted windows at offsets c in {0, +/- L, +/- 2L, ..., +/- 1/4}.
  - v9's Fejer-kernel test functional on the lifted M[k,k].
  - v5/v7's even-symmetry constraints (E)/(F)/(G).
  - The PSD lift on (p, Q) and the standard mass / nonneg constraints.

The relaxation bounds C^{even}_{6.2}.  The new windows do not change
the validity of the Fejer test or the (G) constraint, since they are
all derived without any approximation.
"""

from __future__ import annotations

from dataclasses import dataclass

import cvxpy as cp
import numpy as np


@dataclass
class V11Result:
    status: str
    Omega: float
    primal: float


class AutocorrLowerBoundV11:
    def __init__(self, N: int = 16, K: int = 32) -> None:
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
        self.M = cp.Variable((K + 1, K + 1), symmetric=True, name="M")

        constraints: list = []

        # ---- (P, Q) lift ----
        block = cp.bmat([
            [np.array([[1.0]]), cp.reshape(self.p, (1, self.dim), order="C")],
            [cp.reshape(self.p, (self.dim, 1), order="C"), self.Q],
        ])
        constraints.append(block >> 0)
        constraints.append(cp.sum(self.p) == 1)
        constraints.append(self.Q >= 0)
        constraints.append(cp.sum(self.Q, axis=1) == self.p)
        constraints.append(self.a[0] == 1)

        # even symmetry
        for i in range(self.dim // 2):
            constraints.append(self.p[i] == self.p[self.dim - 1 - i])
        for i in range(self.dim):
            for l in range(self.dim):
                if i + l < self.dim - 1:
                    i2 = self.dim - 1 - i
                    l2 = self.dim - 1 - l
                    if (i2, l2) != (i, l):
                        constraints.append(self.Q[i, l] == self.Q[i2, l2])
        for i in range(self.dim):
            constraints.append(self.Q[i, self.dim - 1 - i] == self.Q[i, i])

        # ---- M PSD lift on (1, a_1, ..., a_K) ----
        constraints.append(self.M >> 0)
        constraints.append(self.M[0, 0] == 1)
        for k in range(1, K + 1):
            constraints.append(self.M[0, k] == self.a[k])
            constraints.append(self.M[k, 0] == self.a[k])

        # f-side linear bounds on a_k from p
        cos_min = np.zeros((K + 1, self.dim))
        cos_max = np.zeros((K + 1, self.dim))
        for k in range(1, K + 1):
            for j in range(self.dim):
                a_left = -0.25 + j * self.L
                a_right = a_left + self.L
                xs = np.linspace(a_left, a_right, 401)
                cv = np.cos(2.0 * np.pi * k * xs)
                cos_min[k, j] = cv.min()
                cos_max[k, j] = cv.max()
        for k in range(1, K + 1):
            constraints.append(self.a[k] >= cos_min[k] @ self.p)
            constraints.append(self.a[k] <= cos_max[k] @ self.p)

        # ---- Fejer kernel constraint ----
        weights = np.array([1.0 - k / (K + 1) for k in range(1, K + 1)])
        m_diags = cp.hstack([self.M[k, k] for k in range(1, K + 1)])
        constraints.append(self.Omega >= 1.0 + 2.0 * (weights @ m_diags))

        # ---- v7 (W) family of windows: centered + shifted ----
        widths = list(np.linspace(2.0 * self.L, 0.5, max(8, self.dim)))  # half-widths
        offsets = [0.0] + [c * self.L for c in range(1, self.dim // 2 + 1)] + [-c * self.L for c in range(1, self.dim // 2 + 1)]
        for c in offsets:
            for h in widths:
                a_w = c - h
                b_w = c + h
                if a_w < -0.5 or b_w > 0.5:
                    continue
                # cells fully inside: 2 + (a_w + 0.5)/L  <=  j+k  <=  (b_w + 0.5)/L
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

    def solve(self, solver: str = "MOSEK", verbose: bool = False, **kwargs) -> V11Result:
        val = self.problem.solve(solver=solver, verbose=verbose, **kwargs)
        return V11Result(
            status=self.problem.status,
            Omega=float(self.Omega.value) if self.Omega.value is not None else float("nan"),
            primal=float(val) if val is not None else float("nan"),
        )


if __name__ == "__main__":
    for N in (8, 16, 24):
        for K in (16, 32):
            try:
                prob = AutocorrLowerBoundV11(N=N, K=K)
                out = prob.solve(solver="MOSEK", verbose=False)
                print(f"N={N:2d}, K={K:2d}: Omega={out.Omega:.6f}")
            except Exception as e:
                print(f"N={N}, K={K}: failed: {e}")
