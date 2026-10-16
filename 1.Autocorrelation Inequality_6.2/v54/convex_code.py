"""v54 = v51 + Putinar localizer h_2(x) = cos(2 pi x) - cos^2(2 pi x).

This localizer is nonneg on [-1/4, 1/4] and negative on (1/4, 3/4),
same support as the v51 localizer, but algebraically distinct.

Fourier coefficients of h_2:
    h_2(x) = -1/2 + cos(2 pi x) - (1/2) cos(4 pi x)
    h_2_hat(0)  = -1/2
    h_2_hat(±1) = +1/2
    h_2_hat(±2) = -1/4
"""

from __future__ import annotations

from dataclasses import dataclass

import cvxpy as cp
import numpy as np


@dataclass
class V54Result:
    status: str
    Omega: float
    primal: float


def _sign(k: int) -> float:
    if k > 0:
        return 1.0
    if k < 0:
        return -1.0
    return 0.0


class AutocorrLowerBoundV54:
    def __init__(self, N: int = 8, K: int = 8) -> None:
        self.N = N
        self.K = K
        self.K_ext = K + 2
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

        # Cell bounds on a_k, b_k for k=1..K+2
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

        for k in range(1, K + 1):
            sum_diag = self.M[k, k] + self.M[K + k, K + k]
            constraints.append(self.v[k - 1] >= cp.square(sum_diag))
        weights = np.array([1.0 - k / (K + 1) for k in range(1, K + 1)])
        constraints.append(self.Omega >= 1.0 + 2.0 * (weights @ self.v))

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

        # Plain Bochner on (a, b) at size K+3
        size_T = self.K_ext + 1
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

        # v51 localized Bochner h_1 = cos(2 pi x), size K+1
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

        # NEW v54: localized Bochner h_2 = cos - cos^2, size K+1
        # h_2_hat(0) = -1/2, h_2_hat(±1) = +1/2, h_2_hat(±2) = -1/4
        # nu_2_hat(k) = -1/2 f_hat(k) + (1/2)(f_hat(k-1) + f_hat(k+1))
        #               -(1/4)(f_hat(k-2) + f_hat(k+2))
        def nu2_re(k: int):
            return (
                -0.5 * self.a[abs(k)]
                + 0.5 * (self.a[abs(k - 1)] + self.a[abs(k + 1)])
                - 0.25 * (self.a[abs(k - 2)] + self.a[abs(k + 2)])
            )

        def nu2_im(k: int):
            # Im(f_hat(m)) = -sign(m) * b[|m|]
            t0 = -_sign(k) * self.b[abs(k)]
            t1 = -_sign(k - 1) * self.b[abs(k - 1)] + -_sign(k + 1) * self.b[abs(k + 1)]
            t2 = -_sign(k - 2) * self.b[abs(k - 2)] + -_sign(k + 2) * self.b[abs(k + 2)]
            return -0.5 * t0 + 0.5 * t1 - 0.25 * t2

        Tn2_R_rows = []
        Tn2_I_rows = []
        for i in range(K + 1):
            row_R = []
            row_I = []
            for j in range(K + 1):
                k = i - j
                row_R.append(nu2_re(k))
                row_I.append(nu2_im(k))
            Tn2_R_rows.append(cp.hstack(row_R))
            Tn2_I_rows.append(cp.hstack(row_I))
        Tn2_R = cp.vstack(Tn2_R_rows)
        Tn2_I = cp.vstack(Tn2_I_rows)
        big_nu2 = cp.bmat([[Tn2_R, -Tn2_I], [Tn2_I, Tn2_R]])
        constraints.append(big_nu2 >> 0)

        self.constraints = constraints
        self.problem = cp.Problem(cp.Minimize(self.Omega), constraints)

    def solve(self, solver: str = "MOSEK", verbose: bool = False, **kwargs) -> V54Result:
        val = self.problem.solve(solver=solver, verbose=verbose, **kwargs)
        return V54Result(
            status=self.problem.status,
            Omega=float(self.Omega.value) if self.Omega.value is not None else float("nan"),
            primal=float(val) if val is not None else float("nan"),
        )


if __name__ == "__main__":
    import time

    print("v54 (v51 + cos-cos^2 localizer) sanity sweep:")
    for N, K in [(8, 8), (16, 16), (16, 32), (24, 24), (32, 32), (50, 32)]:
        t = time.time()
        try:
            prob = AutocorrLowerBoundV54(N=N, K=K)
            out = prob.solve(solver="MOSEK", verbose=False)
            print(f"N={N:3d}, K={K:3d}: Omega={out.Omega:.6f}  ({time.time()-t:.1f}s)")
        except Exception as e:
            print(f"N={N}, K={K}: failed: {e}")
