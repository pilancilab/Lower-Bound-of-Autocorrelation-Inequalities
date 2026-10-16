"""v50 = v45 + extended (a, b) to 2K + Gram matrix G >> 0 in real trig basis.

G[i, j] = integral f(x) g_i(x) g_j(x) dx
where g_0 = 1, g_k = cos(2pi k x), g_{K+k} = sin(2pi k x) for k=1..K.

Via trig product-to-sum identities, the entries of G are LINEAR in
(a_0..2K, b_0..2K):
    G[0, 0] = 1
    G[0, i] = a_i, G[0, K+i] = b_i (for i=1..K)
    G[i, j] = (a_{|i-j|} + a_{i+j}) / 2       (cos-cos)
    G[K+i, K+j] = (a_{|i-j|} - a_{i+j}) / 2   (sin-sin)
    G[i, K+j] = (b_{i+j} + sgn(j-i) * b_{|j-i|}) / 2  (cos-sin)
The constraint is G >> 0.
"""

from __future__ import annotations

from dataclasses import dataclass

import cvxpy as cp
import numpy as np


@dataclass
class V50Result:
    status: str
    Omega: float
    primal: float


class AutocorrLowerBoundV50:
    def __init__(self, N: int = 8, K: int = 8) -> None:
        self.N = N
        self.K = K
        self.K2 = 2 * K
        self.dim_p = 2 * N
        self.L_f = 1.0 / (4 * N)

        self.Omega = cp.Variable(nonneg=True, name="Omega")
        self.p = cp.Variable(self.dim_p, nonneg=True, name="p")
        self.a = cp.Variable(self.K2 + 1, name="a")
        self.b = cp.Variable(self.K2 + 1, name="b")
        self.M = cp.Variable((2 * K + 1, 2 * K + 1), symmetric=True, name="M")
        self.v = cp.Variable(K, nonneg=True, name="v")
        self.Q = cp.Variable((K + 1, K + 1), hermitian=True, name="Q")
        # NEW v50: G in real trig basis (size 2K+1)
        self.G = cp.Variable((2 * K + 1, 2 * K + 1), symmetric=True, name="G")

        constraints: list = []

        constraints.append(cp.sum(self.p) == 1)
        constraints.append(self.a[0] == 1)
        constraints.append(self.b[0] == 0)

        # M lift
        constraints.append(self.M >> 0)
        constraints.append(self.M[0, 0] == 1)
        for k in range(1, K + 1):
            constraints.append(self.M[0, k] == self.a[k])
            constraints.append(self.M[k, 0] == self.a[k])
            constraints.append(self.M[0, K + k] == self.b[k])
            constraints.append(self.M[K + k, 0] == self.b[k])

        # Cell bounds for k=1..2K
        cos_min = np.zeros((self.K2 + 1, self.dim_p))
        cos_max = np.zeros((self.K2 + 1, self.dim_p))
        sin_min = np.zeros((self.K2 + 1, self.dim_p))
        sin_max = np.zeros((self.K2 + 1, self.dim_p))
        for k in range(1, self.K2 + 1):
            for j in range(self.dim_p):
                a_l = -0.25 + j * self.L_f
                a_r = a_l + self.L_f
                xs = np.linspace(a_l, a_r, 401)
                cv = np.cos(2 * np.pi * k * xs)
                sv = np.sin(2 * np.pi * k * xs)
                cos_min[k, j] = cv.min()
                cos_max[k, j] = cv.max()
                sin_min[k, j] = sv.min()
                sin_max[k, j] = sv.max()
        for k in range(1, self.K2 + 1):
            constraints.append(self.a[k] >= cos_min[k] @ self.p)
            constraints.append(self.a[k] <= cos_max[k] @ self.p)
            constraints.append(self.b[k] >= sin_min[k] @ self.p)
            constraints.append(self.b[k] <= sin_max[k] @ self.p)

        # v14 squared (redundant)
        for k in range(1, K + 1):
            sum_diag = self.M[k, k] + self.M[K + k, K + k]
            constraints.append(self.v[k - 1] >= cp.square(sum_diag))
        weights = np.array([1.0 - k / (K + 1) for k in range(1, K + 1)])
        constraints.append(self.Omega >= 1.0 + 2.0 * (weights @ self.v))

        # Fejer-Riesz
        constraints.append(self.Q >> 0)
        constraints.append(cp.real(cp.trace(self.Q)) == self.Omega - 1)
        for k in range(1, K + 1):
            anti_diag_sum = sum(self.Q[i, i + k] for i in range(K + 1 - k))
            constraints.append(
                cp.real(anti_diag_sum)
                == -weights[k - 1] * (self.M[k, k] - self.M[K + k, K + k])
            )
            constraints.append(
                cp.imag(anti_diag_sum) == 2.0 * weights[k - 1] * self.M[k, K + k]
            )

        # Bochner-Toeplitz at size 2K+1 (v46 extension)
        T_R_rows = []
        T_I_rows = []
        for i in range(self.K2 + 1):
            row_R = []
            row_I = []
            for j in range(self.K2 + 1):
                d = abs(i - j)
                if d == 0:
                    row_R.append(self.a[0])
                    row_I.append(cp.Constant(0.0))
                else:
                    row_R.append(self.a[d])
                    if i > j:
                        row_I.append(-self.b[d])
                    else:
                        row_I.append(self.b[d])
            T_R_rows.append(cp.hstack(row_R))
            T_I_rows.append(cp.hstack(row_I))
        T_R = cp.vstack(T_R_rows)
        T_I = cp.vstack(T_I_rows)
        big = cp.bmat([[T_R, -T_I], [T_I, T_R]])
        constraints.append(big >> 0)

        # NEW v50: Gram matrix G in the real trig basis
        # G has size (2K+1) x (2K+1), symmetric, PSD
        constraints.append(self.G >> 0)
        # G[0, 0] = 1
        constraints.append(self.G[0, 0] == 1)
        # G[0, k] = a_k for k=1..K; G[0, K+k] = b_k for k=1..K
        for k in range(1, K + 1):
            constraints.append(self.G[0, k] == self.a[k])
            constraints.append(self.G[k, 0] == self.a[k])
            constraints.append(self.G[0, K + k] == self.b[k])
            constraints.append(self.G[K + k, 0] == self.b[k])
        # cos-cos block: G[i,j] = (a_{|i-j|} + a_{i+j}) / 2
        for i in range(1, K + 1):
            for j in range(1, K + 1):
                d = abs(i - j)
                s = i + j  # <= 2K
                constraints.append(
                    self.G[i, j] == (self.a[d] + self.a[s]) / 2.0
                )
        # sin-sin block: G[K+i, K+j] = (a_{|i-j|} - a_{i+j}) / 2
        for i in range(1, K + 1):
            for j in range(1, K + 1):
                d = abs(i - j)
                s = i + j
                constraints.append(
                    self.G[K + i, K + j] == (self.a[d] - self.a[s]) / 2.0
                )
        # cos-sin block: G[i, K+j] = (b_{i+j} + sign(j-i) * b_{|j-i|}) / 2
        for i in range(1, K + 1):
            for j in range(1, K + 1):
                d = abs(j - i)
                s = i + j
                if j > i:
                    sign = 1.0
                elif j < i:
                    sign = -1.0
                else:
                    sign = 0.0
                expr = (self.b[s] + sign * self.b[d]) / 2.0
                constraints.append(self.G[i, K + j] == expr)
                constraints.append(self.G[K + j, i] == expr)

        self.constraints = constraints
        self.problem = cp.Problem(cp.Minimize(self.Omega), constraints)

    def solve(self, solver: str = "MOSEK", verbose: bool = False, **kwargs) -> V50Result:
        val = self.problem.solve(solver=solver, verbose=verbose, **kwargs)
        return V50Result(
            status=self.problem.status,
            Omega=float(self.Omega.value) if self.Omega.value is not None else float("nan"),
            primal=float(val) if val is not None else float("nan"),
        )


if __name__ == "__main__":
    import time

    print("v50 (v45 + Gram-matrix G) sanity sweep:")
    for N, K in [(8, 8), (16, 16), (16, 32), (24, 24), (32, 32), (50, 32)]:
        t = time.time()
        try:
            prob = AutocorrLowerBoundV50(N=N, K=K)
            out = prob.solve(solver="MOSEK", verbose=False)
            print(f"N={N:3d}, K={K:3d}: Omega={out.Omega:.6f}  ({time.time()-t:.1f}s)")
        except Exception as e:
            print(f"N={N}, K={K}: failed: {e}")
