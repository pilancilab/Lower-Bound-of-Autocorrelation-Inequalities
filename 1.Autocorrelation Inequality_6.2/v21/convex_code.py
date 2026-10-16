"""v21: v14 with SQUARED Fejer kernel as test functional.

Per v19 verdict's recommendation: use a different positive trig
polynomial.  The SQUARED Fejer kernel  P(t) = F_K(t)^2 / Z_K
(normalized to integral 1) has Fourier coefficients that are the
DISCRETE AUTOCORRELATION of the Fejer Fourier coefficients,
specifically:

   P_hat(m)  =  (1/Z_K) * sum_{k} F_K_hat(k) * F_K_hat(k+m)
              =  (1/Z_K) * sum_{k} (1 - |k|/(K+1)) * (1 - |k+m|/(K+1)),

valid for |m| <= 2K.  Z_K = sum_{k} (1 - |k|/(K+1))^2 (the L^2 norm
of F_K, also = P_hat at m = 0 before normalization).

P is nonneg (it's the square of a real function), and integrates to
1 by construction.  So we can use P in place of F_K in the v14
constraint:

   1 + 2 sum_{m=1}^{2K} P_hat(m) * |F_hat(m)|^2  <=  Omega.

The hope is that the squared Fejer's broader spectral support
(|m| <= 2K instead of |m| <= K) and different weighting captures
more information about |F_hat(m)|^2.
"""

from __future__ import annotations

from dataclasses import dataclass

import cvxpy as cp
import numpy as np


@dataclass
class V21Result:
    status: str
    Omega: float
    primal: float


def squared_fejer_coeffs(K: int):
    """Return P_hat(m) for m = 0, 1, ..., 2K (the rest by symmetry)."""
    fejer = np.array([1.0 - abs(k) / (K + 1) for k in range(-K, K + 1)])
    coeffs = np.convolve(fejer, fejer, mode="full")  # length 2*(2K+1)-1 = 4K+1
    # autocorrelation of fejer = convolution of fejer with itself reversed.
    # For the Fejer kernel with real symmetric coeffs, autocorrelation = convolution.
    Z_K = coeffs[len(coeffs) // 2]  # central value
    coeffs = coeffs / Z_K
    # Extract coeffs at positions m = 0, 1, 2, ..., 2K
    center = len(coeffs) // 2
    return np.array([coeffs[center + m] for m in range(0, 2 * K + 1)])


class AutocorrLowerBoundV21:
    def __init__(self, N: int = 8, K: int = 8) -> None:
        if N < 1 or K < 1:
            raise ValueError("N, K must be positive")
        self.N = N
        self.K = K
        self.K_max = 2 * K  # squared Fejer has degree 2K
        self.dim_p = 2 * N
        self.L_f = 1.0 / (4 * N)

        self.Omega = cp.Variable(nonneg=True, name="Omega")
        self.p = cp.Variable(self.dim_p, nonneg=True, name="p")
        self.a = cp.Variable(self.K_max + 1, name="a")
        self.b = cp.Variable(self.K_max + 1, name="b")
        self.M = cp.Variable((2 * self.K_max + 1, 2 * self.K_max + 1), symmetric=True, name="M")
        self.v = cp.Variable(self.K_max, nonneg=True, name="v")

        constraints: list = []

        constraints.append(cp.sum(self.p) == 1)
        constraints.append(self.a[0] == 1)
        constraints.append(self.b[0] == 0)

        constraints.append(self.M >> 0)
        constraints.append(self.M[0, 0] == 1)
        for k in range(1, self.K_max + 1):
            constraints.append(self.M[0, k] == self.a[k])
            constraints.append(self.M[k, 0] == self.a[k])
            constraints.append(self.M[0, self.K_max + k] == self.b[k])
            constraints.append(self.M[self.K_max + k, 0] == self.b[k])

        cos_min = np.zeros((self.K_max + 1, self.dim_p))
        cos_max = np.zeros((self.K_max + 1, self.dim_p))
        sin_min = np.zeros((self.K_max + 1, self.dim_p))
        sin_max = np.zeros((self.K_max + 1, self.dim_p))
        for k in range(1, self.K_max + 1):
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
        for k in range(1, self.K_max + 1):
            constraints.append(self.a[k] >= cos_min[k] @ self.p)
            constraints.append(self.a[k] <= cos_max[k] @ self.p)
            constraints.append(self.b[k] >= sin_min[k] @ self.p)
            constraints.append(self.b[k] <= sin_max[k] @ self.p)

        for k in range(1, self.K_max + 1):
            sum_diag = self.M[k, k] + self.M[self.K_max + k, self.K_max + k]
            constraints.append(self.v[k - 1] >= cp.square(sum_diag))

        # Squared Fejer kernel coefficients
        sq_coeffs = squared_fejer_coeffs(K)  # length 2K+1, indices 0..2K
        # constraint: Omega >= 1 + 2 sum_{m=1}^{2K} sq_coeffs[m] * v_{m-1}
        weights_sq = sq_coeffs[1:]
        v_slice = cp.hstack([self.v[m - 1] for m in range(1, self.K_max + 1)])
        constraints.append(self.Omega >= 1.0 + 2.0 * (weights_sq @ v_slice))

        self.constraints = constraints
        self.problem = cp.Problem(cp.Minimize(self.Omega), constraints)

    def solve(self, solver: str = "MOSEK", verbose: bool = False, **kwargs) -> V21Result:
        val = self.problem.solve(solver=solver, verbose=verbose, **kwargs)
        return V21Result(
            status=self.problem.status,
            Omega=float(self.Omega.value) if self.Omega.value is not None else float("nan"),
            primal=float(val) if val is not None else float("nan"),
        )


if __name__ == "__main__":
    for N in (8, 16):
        for K in (4, 8, 16):
            try:
                prob = AutocorrLowerBoundV21(N=N, K=K)
                out = prob.solve(solver="MOSEK", verbose=False)
                print(f"N={N:2d}, K={K:2d} (eff degree {2*K}): Omega={out.Omega:.6f}")
            except Exception as e:
                print(f"N={N}, K={K}: failed: {e}")
