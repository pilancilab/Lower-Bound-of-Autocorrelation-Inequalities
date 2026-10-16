"""v51_verified = v51 with ANALYTICAL cell bounds (no grid error).

The original v51 used a 401-point grid to estimate min/max of
cos(2 pi k x) and sin(2 pi k x) over each cell [a_l, a_r]. The grid
min is >= the true min (because min over a sample >= min over a
superset), which makes the cell-bound constraint a_k >= cos_min · p
STRICTER than mathematically correct, possibly excluding admissible f.

This version computes min/max analytically:
- cos(2 pi k x) has critical points at x = m/(2k), with value (-1)^m.
- sin(2 pi k x) has critical points at x = (2m+1)/(4k), with value (-1)^m.

For a cell [a_l, a_r], the true extremum is at either an endpoint
or at any critical point lying strictly inside the cell.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import ceil, floor

import cvxpy as cp
import numpy as np


@dataclass
class V51vResult:
    status: str
    Omega: float
    primal: float


def _sign(k: int) -> float:
    if k > 0:
        return 1.0
    if k < 0:
        return -1.0
    return 0.0


def cos_extrema_on_cell(k: int, a_l: float, a_r: float) -> tuple[float, float]:
    """Compute exact min and max of cos(2 pi k x) on [a_l, a_r].

    cos(2 pi k x) has critical points at x = m/(2k) for integer m,
    with value (-1)^m at that point. (This is because d/dx cos(2pi k x) =
    -2pi k sin(2pi k x) = 0 when 2 pi k x = m pi, i.e., x = m/(2k).)

    The extremum on [a_l, a_r] is at either endpoint or at a critical
    point strictly inside the cell.
    """
    cand_vals = [np.cos(2.0 * np.pi * k * a_l), np.cos(2.0 * np.pi * k * a_r)]
    # Critical points: x = m/(2k) for integer m.
    # Find all m with m/(2k) in (a_l, a_r).
    m_lo = int(ceil(2.0 * k * a_l + 1e-12))  # strict > a_l
    m_hi = int(floor(2.0 * k * a_r - 1e-12))  # strict < a_r
    for m in range(m_lo, m_hi + 1):
        x_c = m / (2.0 * k)
        if a_l < x_c < a_r:
            cand_vals.append((-1.0) ** m)  # cos(pi m)
    return min(cand_vals), max(cand_vals)


def sin_extrema_on_cell(k: int, a_l: float, a_r: float) -> tuple[float, float]:
    """Compute exact min and max of sin(2 pi k x) on [a_l, a_r].

    sin(2 pi k x) has critical points at x where cos(2 pi k x) = 0,
    i.e., 2 pi k x = (m + 1/2) pi for integer m, i.e., x = (2m+1)/(4k).
    At these points, sin(2 pi k x) = sin((m + 1/2) pi) = (-1)^m.
    """
    cand_vals = [np.sin(2.0 * np.pi * k * a_l), np.sin(2.0 * np.pi * k * a_r)]
    # Critical points: x = (2m+1)/(4k) for integer m.
    # Find all m with (2m+1)/(4k) in (a_l, a_r).
    # (2m+1)/(4k) > a_l => m > 2 k a_l - 1/2
    # (2m+1)/(4k) < a_r => m < 2 k a_r - 1/2
    m_lo = int(ceil(2.0 * k * a_l - 0.5 + 1e-12))
    m_hi = int(floor(2.0 * k * a_r - 0.5 - 1e-12))
    for m in range(m_lo, m_hi + 1):
        x_c = (2 * m + 1) / (4.0 * k)
        if a_l < x_c < a_r:
            cand_vals.append((-1.0) ** m)
    return min(cand_vals), max(cand_vals)


class AutocorrLowerBoundV51v:
    def __init__(self, N: int = 8, K: int = 8) -> None:
        self.N = N
        self.K = K
        self.K_ext = K + 1
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

        # ANALYTICAL cell bounds on a_k, b_k for k=1..K+1
        cos_min = np.zeros((self.K_ext + 1, self.dim_p))
        cos_max = np.zeros((self.K_ext + 1, self.dim_p))
        sin_min = np.zeros((self.K_ext + 1, self.dim_p))
        sin_max = np.zeros((self.K_ext + 1, self.dim_p))
        for k in range(1, self.K_ext + 1):
            for j in range(self.dim_p):
                a_l = -0.25 + j * self.L_f
                a_r = a_l + self.L_f
                cmn, cmx = cos_extrema_on_cell(k, a_l, a_r)
                smn, smx = sin_extrema_on_cell(k, a_l, a_r)
                cos_min[k, j] = cmn
                cos_max[k, j] = cmx
                sin_min[k, j] = smn
                sin_max[k, j] = smx
        for k in range(1, self.K_ext + 1):
            constraints.append(self.a[k] >= cos_min[k] @ self.p)
            constraints.append(self.a[k] <= cos_max[k] @ self.p)
            constraints.append(self.b[k] >= sin_min[k] @ self.p)
            constraints.append(self.b[k] <= sin_max[k] @ self.p)

        # v14 squared (loose, redundant)
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

        # Plain Bochner-Toeplitz size K+1
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

        # Localized Bochner with h(x) = cos(2 pi x)
        def nu_re(k: int):
            return (self.a[abs(k - 1)] + self.a[abs(k + 1)]) / 2.0

        def nu_im(k: int):
            t1 = -_sign(k - 1) * self.b[abs(k - 1)]
            t2 = -_sign(k + 1) * self.b[abs(k + 1)]
            return (t1 + t2) / 2.0

        Tn_R_rows = []
        Tn_I_rows = []
        for i in range(K + 1):
            row_R = []
            row_I = []
            for j in range(K + 1):
                k = i - j
                row_R.append(nu_re(k))
                row_I.append(nu_im(k))
            Tn_R_rows.append(cp.hstack(row_R))
            Tn_I_rows.append(cp.hstack(row_I))
        Tn_R = cp.vstack(Tn_R_rows)
        Tn_I = cp.vstack(Tn_I_rows)
        big_nu = cp.bmat([[Tn_R, -Tn_I], [Tn_I, Tn_R]])
        constraints.append(big_nu >> 0)

        self.constraints = constraints
        self.problem = cp.Problem(cp.Minimize(self.Omega), constraints)

    def solve(self, solver: str = "MOSEK", verbose: bool = False, **kwargs) -> V51vResult:
        val = self.problem.solve(solver=solver, verbose=verbose, **kwargs)
        return V51vResult(
            status=self.problem.status,
            Omega=float(self.Omega.value) if self.Omega.value is not None else float("nan"),
            primal=float(val) if val is not None else float("nan"),
        )


if __name__ == "__main__":
    import time

    # First, sanity check the analytical extrema against the grid
    print("Sanity check: analytical vs grid extrema (should agree to ~1e-7):")
    for k in [1, 8, 32, 64, 97]:
        a_l, a_r = -0.25 + 0.5 * 1.0 / 8.0, -0.25 + 0.5 * 2.0 / 8.0  # some test cell
        cmn_a, cmx_a = cos_extrema_on_cell(k, a_l, a_r)
        xs = np.linspace(a_l, a_r, 4001)  # much finer grid
        cv = np.cos(2 * np.pi * k * xs)
        cmn_g, cmx_g = cv.min(), cv.max()
        print(f"  k={k}: analytical=[{cmn_a:.8f}, {cmx_a:.8f}], grid=[{cmn_g:.8f}, {cmx_g:.8f}], diff_min={cmn_a-cmn_g:+.2e}")

    print("\nv51_verified sanity sweep:")
    for N, K in [(8, 8), (16, 16), (16, 32), (32, 32), (50, 32)]:
        t = time.time()
        try:
            prob = AutocorrLowerBoundV51v(N=N, K=K)
            out = prob.solve(solver="MOSEK", verbose=False)
            print(f"N={N:3d}, K={K:3d}: Omega={out.Omega:.6f}  ({time.time()-t:.1f}s)")
        except Exception as e:
            print(f"N={N}, K={K}: failed: {e}")
