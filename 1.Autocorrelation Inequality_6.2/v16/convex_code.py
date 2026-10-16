"""v16: v14 autocorrelation Fejer + v13 Bochner Toeplitz + v11 windows.

Combines all major rigorous constraints discovered so far.

Constraints:
* (p, Q) PSD lift, mass, marginals, Q >= 0.
* (W) family of windowed integral inequalities (centered + shifted).
* (a, b) f-side Fourier coefficients with linear bounds from p.
* M PSD lift on (1, a_1..a_K, b_1..b_K).
* v14 autocorrelation Fejer constraint:
    Omega >= 1 + 2 * sum_{k=1}^K (1 - k/(K+1)) * v_k
    with rotated SOC v_k >= (M[k,k] + M[K+k, K+k])^2.
* v13 Bochner Toeplitz on the f-side Fourier sequence (a_k - i b_k).

All constraints are individually valid for any admissible f
(no even restriction).
"""

from __future__ import annotations

from dataclasses import dataclass

import cvxpy as cp
import numpy as np


@dataclass
class V16Result:
    status: str
    Omega: float
    primal: float


class AutocorrLowerBoundV16:
    def __init__(self, N: int = 16, K: int = 16, with_windows: bool = True) -> None:
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

        # (p, Q) PSD lift
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

        # M PSD lift
        constraints.append(self.M >> 0)
        constraints.append(self.M[0, 0] == 1)
        for k in range(1, K + 1):
            constraints.append(self.M[0, k] == self.a[k])
            constraints.append(self.M[k, 0] == self.a[k])
            constraints.append(self.M[0, K + k] == self.b[k])
            constraints.append(self.M[K + k, 0] == self.b[k])

        # f-side bounds
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

        # v14 autocorrelation Fejer
        for k in range(1, K + 1):
            sum_diag = self.M[k, k] + self.M[K + k, K + k]
            constraints.append(self.v[k - 1] >= cp.square(sum_diag))
        weights = np.array([1.0 - k / (K + 1) for k in range(1, K + 1)])
        constraints.append(self.Omega >= 1.0 + 2.0 * (weights @ self.v))

        # v13 Bochner Toeplitz on (a, b)
        T_R_rows = []
        T_I_rows = []
        for i in range(K + 1):
            row_R = []
            row_I = []
            for j in range(K + 1):
                d = abs(i - j)
                if d == 0:
                    row_R.append(self.a[0])
                    row_I.append(cp.Constant(0.0))
                else:
                    row_R.append(self.a[d])
                    sign = -1.0 if i > j else 1.0
                    row_I.append(sign * self.b[d])
            T_R_rows.append(cp.hstack(row_R))
            T_I_rows.append(cp.hstack(row_I))
        T_R = cp.vstack(T_R_rows)
        T_I = cp.vstack(T_I_rows)
        big = cp.bmat([[T_R, -T_I], [T_I, T_R]])
        constraints.append(big >> 0)

        # v11 (W) windowed family (optional, heavy)
        if with_windows:
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

    def solve(self, solver: str = "MOSEK", verbose: bool = False, **kwargs) -> V16Result:
        val = self.problem.solve(solver=solver, verbose=verbose, **kwargs)
        return V16Result(
            status=self.problem.status,
            Omega=float(self.Omega.value) if self.Omega.value is not None else float("nan"),
            primal=float(val) if val is not None else float("nan"),
        )


if __name__ == "__main__":
    for N, K, w in [(8, 8, True), (16, 16, True), (24, 24, False), (32, 32, False)]:
        try:
            prob = AutocorrLowerBoundV16(N=N, K=K, with_windows=w)
            out = prob.solve(solver="MOSEK", verbose=False)
            print(f"N={N:2d}, K={K:2d}, windows={w}: Omega={out.Omega:.6f}")
        except Exception as e:
            print(f"N={N}, K={K}: failed: {e}")
