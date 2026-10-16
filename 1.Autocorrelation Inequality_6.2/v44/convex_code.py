"""v44: v30 with Fejer-Riesz exact trig-poly nonnegativity.

Replaces v30's n_tests grid sampling of sigma_K^{f*f}(t*) <= Omega
with a single Fejer-Riesz SDP constraint enforcing
   sup_t sigma_K^{f*f}(t) <= Omega
exactly.

The constraint is:
    Omega - sigma_K^{f*f}(t) >= 0  for all t
which for a degree-K trig polynomial is equivalent to existence of
a Hermitian PSD matrix Q in C^{(K+1)x(K+1)} such that
    sum_{i=0}^{K-k} Q[i, i+k] = c_k for k = 0, 1, ..., K
where
    c_0 = Omega - 1
    c_k = -w_k * (M[k,k] - M[K+k,K+k]) - 2i * w_k * M[k, K+k]
    w_k = 1 - k/(K+1).

The v14 autocorrelation Fejer constraint and all other v30
constraints are inherited unchanged.
"""

from __future__ import annotations

from dataclasses import dataclass

import cvxpy as cp
import numpy as np


@dataclass
class V44Result:
    status: str
    Omega: float
    primal: float


class AutocorrLowerBoundV44:
    def __init__(self, N: int = 8, K: int = 8) -> None:
        self.N = N
        self.K = K
        self.dim_p = 2 * N
        self.L_f = 1.0 / (4 * N)

        self.Omega = cp.Variable(nonneg=True, name="Omega")
        self.p = cp.Variable(self.dim_p, nonneg=True, name="p")
        self.a = cp.Variable(K + 1, name="a")
        self.b = cp.Variable(K + 1, name="b")
        self.M = cp.Variable((2 * K + 1, 2 * K + 1), symmetric=True, name="M")
        self.v = cp.Variable(K, nonneg=True, name="v")
        # Fejer-Riesz Hermitian PSD matrix
        self.Q = cp.Variable((K + 1, K + 1), hermitian=True, name="Q")

        constraints: list = []

        # Mass normalization, a_0 = 1, b_0 = 0
        constraints.append(cp.sum(self.p) == 1)
        constraints.append(self.a[0] == 1)
        constraints.append(self.b[0] == 0)

        # Moment lift: M PSD with first-row equalities
        constraints.append(self.M >> 0)
        constraints.append(self.M[0, 0] == 1)
        for k in range(1, K + 1):
            constraints.append(self.M[0, k] == self.a[k])
            constraints.append(self.M[k, 0] == self.a[k])
            constraints.append(self.M[0, K + k] == self.b[k])
            constraints.append(self.M[K + k, 0] == self.b[k])

        # Cell-bounds on a_k, b_k from p
        cos_min = np.zeros((K + 1, self.dim_p))
        cos_max = np.zeros((K + 1, self.dim_p))
        sin_min = np.zeros((K + 1, self.dim_p))
        sin_max = np.zeros((K + 1, self.dim_p))
        for k in range(1, K + 1):
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
        for k in range(1, K + 1):
            constraints.append(self.a[k] >= cos_min[k] @ self.p)
            constraints.append(self.a[k] <= cos_max[k] @ self.p)
            constraints.append(self.b[k] >= sin_min[k] @ self.p)
            constraints.append(self.b[k] <= sin_max[k] @ self.p)

        # v14 autocorrelation Fejer (kept as redundant safety net)
        for k in range(1, K + 1):
            sum_diag = self.M[k, k] + self.M[K + k, K + k]
            constraints.append(self.v[k - 1] >= cp.square(sum_diag))
        weights = np.array([1.0 - k / (K + 1) for k in range(1, K + 1)])
        constraints.append(self.Omega >= 1.0 + 2.0 * (weights @ self.v))

        # NEW v44: Fejer-Riesz exact constraint
        # Q is Hermitian PSD of size (K+1)x(K+1)
        constraints.append(self.Q >> 0)
        # c_0 = Omega - 1 = sum_i Q[i,i] = trace(Q)
        # trace of a hermitian matrix is real automatically
        constraints.append(cp.real(cp.trace(self.Q)) == self.Omega - 1)
        # c_k = sum_{i=0}^{K-k} Q[i, i+k] for k = 1, ..., K
        for k in range(1, K + 1):
            # k-th anti-diagonal sum (offset k above main diag)
            anti_diag_sum = sum(self.Q[i, i + k] for i in range(K + 1 - k))
            # Real part: -w_k * (M[k,k] - M[K+k, K+k])
            constraints.append(
                cp.real(anti_diag_sum)
                == -weights[k - 1] * (self.M[k, k] - self.M[K + k, K + k])
            )
            # Imag part: -2 * w_k * M[k, K+k]
            #
            # Sign derivation:
            #   c_k = -w_k * (a_k - i b_k)^2
            #       = -w_k * (a_k^2 - b_k^2 - 2i a_k b_k)
            #       = -w_k * (a_k^2 - b_k^2) + 2i w_k a_k b_k
            #
            # So Im(c_k) = +2 w_k * a_k b_k = +2 w_k * M[k, K+k].
            constraints.append(
                cp.imag(anti_diag_sum) == 2.0 * weights[k - 1] * self.M[k, K + k]
            )

        self.constraints = constraints
        self.problem = cp.Problem(cp.Minimize(self.Omega), constraints)

    def solve(self, solver: str = "MOSEK", verbose: bool = False, **kwargs) -> V44Result:
        val = self.problem.solve(solver=solver, verbose=verbose, **kwargs)
        return V44Result(
            status=self.problem.status,
            Omega=float(self.Omega.value) if self.Omega.value is not None else float("nan"),
            primal=float(val) if val is not None else float("nan"),
        )


if __name__ == "__main__":
    import time

    print("v44 (Fejer-Riesz) sanity sweep:")
    for N, K in [(8, 8), (16, 16), (16, 32), (24, 24), (32, 32), (50, 32)]:
        t = time.time()
        try:
            prob = AutocorrLowerBoundV44(N=N, K=K)
            out = prob.solve(solver="MOSEK", verbose=False)
            print(f"N={N:3d}, K={K:3d}: Omega={out.Omega:.6f}  ({time.time()-t:.1f}s)")
        except Exception as e:
            print(f"N={N}, K={K}: failed: {e}")
