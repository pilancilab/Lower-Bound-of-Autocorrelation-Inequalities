"""v53 = v51 + additional Putinar localizer h_3(x) = cos^3(2 pi x)."""

from __future__ import annotations

from dataclasses import dataclass

import cvxpy as cp
import numpy as np


@dataclass
class V53Result:
    status: str
    Omega: float
    primal: float


def _sign(k: int) -> float:
    if k > 0:
        return 1.0
    if k < 0:
        return -1.0
    return 0.0


class AutocorrLowerBoundV53:
    def __init__(self, N: int = 8, K: int = 8) -> None:
        self.N = N
        self.K = K
        # Need a, b up to K+3 for the cubic localizer
        self.K_ext = K + 3
        self.dim_p = 2 * N
        self.L_f = 1.0 / (4 * N)

        self.Omega = cp.Variable(nonneg=True, name="Omega")
        self.p = cp.Variable(self.dim_p, nonneg=True, name="p")
        self.a = cp.Variable(self.K_ext + 1, name="a")
        self.b = cp.Variable(self.K_ext + 1, name="b")
        self.M = cp.Variable((2 * K + 1, 2 * K + 1), symmetric=True, name="M")
        self.v = cp.Variable(K, nonneg=True, name="v")
        self.Q = cp.Variable((K + 1, K + 1), hermitian=True, name="Q")

        constraints: list = []

        constraints.append(cp.sum(self.p) == 1)
        constraints.append(self.a[0] == 1)
        constraints.append(self.b[0] == 0)

        constraints.append(self.M >> 0)
        constraints.append(self.M[0, 0] == 1)
        for k in range(1, K + 1):
            constraints.append(self.M[0, k] == self.a[k])
            constraints.append(self.M[k, 0] == self.a[k])
            constraints.append(self.M[0, K + k] == self.b[k])
            constraints.append(self.M[K + k, 0] == self.b[k])

        # Cell bounds on a_k, b_k for k=1..K+3
        cos_min = np.zeros((self.K_ext + 1, self.dim_p))
        cos_max = np.zeros((self.K_ext + 1, self.dim_p))
        sin_min = np.zeros((self.K_ext + 1, self.dim_p))
        sin_max = np.zeros((self.K_ext + 1, self.dim_p))
        for k in range(1, self.K_ext + 1):
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
        for k in range(1, self.K_ext + 1):
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

        # Fejer-Riesz (degree K)
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

        # Plain Bochner on (a, b) at size K+4 (to use K+3 moments)
        size_T = self.K_ext + 1  # = K+4
        T_R_rows = []
        T_I_rows = []
        for i in range(size_T):
            row_R = []
            row_I = []
            for j in range(size_T):
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

        # Localized Bochner with h_1(x) = cos(2 pi x), size K+1 (from v51)
        def nu1_re(k: int):
            return (self.a[abs(k - 1)] + self.a[abs(k + 1)]) / 2.0

        def nu1_im(k: int):
            t1 = -_sign(k - 1) * self.b[abs(k - 1)]
            t2 = -_sign(k + 1) * self.b[abs(k + 1)]
            return (t1 + t2) / 2.0

        Tn1_R_rows = []
        Tn1_I_rows = []
        for i in range(K + 1):
            row_R = []
            row_I = []
            for j in range(K + 1):
                k = i - j
                row_R.append(nu1_re(k))
                row_I.append(nu1_im(k))
            Tn1_R_rows.append(cp.hstack(row_R))
            Tn1_I_rows.append(cp.hstack(row_I))
        Tn1_R = cp.vstack(Tn1_R_rows)
        Tn1_I = cp.vstack(Tn1_I_rows)
        big_nu1 = cp.bmat([[Tn1_R, -Tn1_I], [Tn1_I, Tn1_R]])
        constraints.append(big_nu1 >> 0)

        # NEW v53: Localized Bochner with h_3(x) = cos^3(2 pi x), size K+1
        # nu_3_hat(k) = (3/8)(a[|k-1|] + a[|k+1|]) + (1/8)(a[|k-3|] + a[|k+3|])
        #               - (3i/8)(sign(k-1) b[|k-1|] + sign(k+1) b[|k+1|])
        #               - (i/8)(sign(k-3) b[|k-3|] + sign(k+3) b[|k+3|])
        def nu3_re(k: int):
            return (
                (3.0 / 8.0) * (self.a[abs(k - 1)] + self.a[abs(k + 1)])
                + (1.0 / 8.0) * (self.a[abs(k - 3)] + self.a[abs(k + 3)])
            )

        def nu3_im(k: int):
            t_main = (
                -_sign(k - 1) * self.b[abs(k - 1)] - _sign(k + 1) * self.b[abs(k + 1)]
            )
            t_high = (
                -_sign(k - 3) * self.b[abs(k - 3)] - _sign(k + 3) * self.b[abs(k + 3)]
            )
            return (3.0 / 8.0) * t_main + (1.0 / 8.0) * t_high

        Tn3_R_rows = []
        Tn3_I_rows = []
        for i in range(K + 1):
            row_R = []
            row_I = []
            for j in range(K + 1):
                k = i - j
                row_R.append(nu3_re(k))
                row_I.append(nu3_im(k))
            Tn3_R_rows.append(cp.hstack(row_R))
            Tn3_I_rows.append(cp.hstack(row_I))
        Tn3_R = cp.vstack(Tn3_R_rows)
        Tn3_I = cp.vstack(Tn3_I_rows)
        big_nu3 = cp.bmat([[Tn3_R, -Tn3_I], [Tn3_I, Tn3_R]])
        constraints.append(big_nu3 >> 0)

        self.constraints = constraints
        self.problem = cp.Problem(cp.Minimize(self.Omega), constraints)

    def solve(self, solver: str = "MOSEK", verbose: bool = False, **kwargs) -> V53Result:
        val = self.problem.solve(solver=solver, verbose=verbose, **kwargs)
        return V53Result(
            status=self.problem.status,
            Omega=float(self.Omega.value) if self.Omega.value is not None else float("nan"),
            primal=float(val) if val is not None else float("nan"),
        )


if __name__ == "__main__":
    import time

    print("v53 (v51 + cos^3 localizer) sanity sweep:")
    for N, K in [(8, 8), (16, 16), (16, 32), (24, 24), (32, 32), (50, 32)]:
        t = time.time()
        try:
            prob = AutocorrLowerBoundV53(N=N, K=K)
            out = prob.solve(solver="MOSEK", verbose=False)
            print(f"N={N:3d}, K={K:3d}: Omega={out.Omega:.6f}  ({time.time()-t:.1f}s)")
        except Exception as e:
            print(f"N={N}, K={K}: failed: {e}")
