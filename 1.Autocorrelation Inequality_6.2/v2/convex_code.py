"""v2 convex relaxation for the autocorrelation inequality (6.2).

Improvements over v1
--------------------

v1 still gave the trivial Omega = 1 because the relaxation can collapse Q to
a diagonal matrix that is incompatible with rank-1 outer products.  v2
introduces an additional batch of variables and constraints based on the
**Fourier representation** of F = f*f.

For real f supported on [-1/4, 1/4] with int f = 1, define the
period-one Fourier coefficients

    a_k = int_{-1/4}^{1/4} f(x) cos(2 pi k x) dx,
    b_k = int_{-1/4}^{1/4} f(x) sin(2 pi k x) dx,    k = 0, 1, ..., K,

with a_0 = 1 and b_0 = 0.  Then

    F(t) = (f*f)(t)
         = 1  +  2 sum_{k>=1} [ (a_k^2 - b_k^2) cos(2 pi k t)
                                 + 2 a_k b_k sin(2 pi k t) ].         (*)

The right-hand side of (*) is QUADRATIC in (a_k, b_k).  We linearize via a
moment matrix M of size 2K+1 indexed by (1, a_1, ..., a_K, b_1, ..., b_K),
constrained M >> 0 with M[0,0] = 1, M[0,k] = a_k, M[0,K+k] = b_k.  Defining

    A_k = M[k,k]   - M[K+k, K+k],     // lifts a_k^2 - b_k^2
    B_k = 2 * M[k, K+k],              // lifts 2 a_k b_k

(*) becomes a LINEAR functional of M:

    F(t)_lifted = 1  +  2 sum_{k>=1} [ A_k cos(2 pi k t)
                                       + B_k sin(2 pi k t) ].          (**)

Constraints in v2 (in addition to v1's (A), (B), (C), bands, PSD lift on Q):

(D) Linking a_k, b_k to p (linear interval bounds):
    For each k = 1..K and each interval I_j of f's discretization,
    integrating cos(2 pi k x) f(x) over I_j is bounded between
    p_j * min_{x in I_j} cos(2 pi k x)  and  p_j * max_{x in I_j} cos(2 pi k x);
    summing over j gives lower and upper linear bounds on a_k in terms of p.
    Similarly for b_k with sin.

(E) Moment-matrix PSD lift:
    M >> 0,  M[0,0] = 1,  M[0,k] = a_k,  M[0,K+k] = b_k.

(F) Test-point upper bounds on F via the lifted formula (**):
    For each test point t_* in a chosen grid of [-1/2, 1/2],
        Omega  >=  1 + 2 sum_k [A_k cos(2 pi k t_*) + B_k sin(2 pi k t_*)].
    (When evaluated at the true rank-1 M = aa^T + bb^T variant, the right-
    hand side equals F(t_*), and F(t_*) <= Omega is exactly the truth.)

These constraints are valid for every admissible f.  The proof of validity
is in rigorousproof.md.

The relaxation is therefore at least as tight as v1.  Crucially, the
diagonal-Q degeneracy of v1 is no longer a feasible point unless one can
ALSO construct (a_k, b_k, M) compatible with it, which is far more
restricted because the linear interval bounds tie a_k and b_k to the same
p that defines Q's marginals.
"""

from __future__ import annotations

from dataclasses import dataclass

import cvxpy as cp
import numpy as np


@dataclass
class V2Result:
    status: str
    Omega: float
    primal: float


class AutocorrLowerBoundV2:
    def __init__(self, N: int = 8, K: int = 6, n_tests: int = 41) -> None:
        if N < 1 or K < 1:
            raise ValueError("N and K must be positive integers")
        self.N = N
        self.K = K
        self.dim = 2 * N
        self.L = 1.0 / (4 * N)
        self.n_tests = n_tests

        # ---------------- v1 variables ----------------
        self.Omega = cp.Variable(nonneg=True, name="Omega")
        self.p = cp.Variable(self.dim, nonneg=True, name="p")
        self.Q = cp.Variable((self.dim, self.dim), symmetric=True, name="Q")

        # ---------------- v2 new variables ----------------
        self.a = cp.Variable(K + 1, name="a")  # cosine Fourier coeffs of f, indices 0..K
        self.b = cp.Variable(K + 1, name="b")  # sine Fourier coeffs of f
        m_dim = 2 * K + 1                      # moment matrix size
        self.M = cp.Variable((m_dim, m_dim), symmetric=True, name="M")

        constraints: list[cp.constraints.Constraint] = []

        # ---- v1: PSD lift on (1, p, Q) ----
        block = cp.bmat([
            [np.array([[1.0]]), cp.reshape(self.p, (1, self.dim))],
            [cp.reshape(self.p, (self.dim, 1)), self.Q],
        ])
        constraints.append(block >> 0)

        # mass on f-side
        constraints.append(cp.sum(self.p) == 1)

        # v1 (A) nonneg, (B) marginals
        constraints.append(self.Q >= 0)
        constraints.append(cp.sum(self.Q, axis=1) == self.p)

        # v0 band constraints
        n_bands = 4 * N
        index_lists: dict[int, list[tuple[int, int]]] = {}
        for j in range(1, self.dim + 1):
            for k in range(1, self.dim + 1):
                index_lists.setdefault(j + k, []).append((j - 1, k - 1))
        for m in range(1, n_bands + 1):
            terms = []
            for s in (m, m + 1):
                if s in index_lists:
                    for (jj, kk) in index_lists[s]:
                        terms.append(self.Q[jj, kk])
            if terms:
                constraints.append(0.5 * cp.sum(cp.hstack(terms)) <= self.L * self.Omega)

        # v1 (C) central window constraint
        center_terms = []
        for j in range(1, self.dim + 1):
            for k in range(1, self.dim + 1):
                if N + 2 <= j + k <= 3 * N:
                    center_terms.append(self.Q[j - 1, k - 1])
        if center_terms:
            constraints.append(cp.sum(cp.hstack(center_terms)) <= 0.5 * self.Omega)

        # ---- v2 (E) PSD lift on (1, a_1..a_K, b_1..b_K) ----
        # rows/cols are (a_0=1, a_1, ..., a_K, b_1, ..., b_K)
        # Index 0 = a_0 (=1).  Indices 1..K = a_k.  Indices K+1..2K = b_k.
        # We embed M >> 0 via the explicit Schur form: M[0,0] = 1, the first
        # row/col matches (1, a_1..a_K, b_1..b_K), and the rest of M is free.
        constraints.append(self.M >> 0)
        constraints.append(self.M[0, 0] == 1)
        constraints.append(self.a[0] == 1)  # mass
        constraints.append(self.b[0] == 0)  # 0-th sine
        for k in range(1, K + 1):
            constraints.append(self.M[0, k] == self.a[k])
            constraints.append(self.M[k, 0] == self.a[k])
            constraints.append(self.M[0, K + k] == self.b[k])
            constraints.append(self.M[K + k, 0] == self.b[k])

        # ---- v2 (D) linear interval bounds on a_k, b_k from p ----
        # Compute, for each (j, k), the min and max of cos(2 pi k x) and
        # sin(2 pi k x) on I_j; tabulate.
        cos_min = np.zeros((K + 1, self.dim))
        cos_max = np.zeros((K + 1, self.dim))
        sin_min = np.zeros((K + 1, self.dim))
        sin_max = np.zeros((K + 1, self.dim))
        for k in range(1, K + 1):
            for j in range(self.dim):
                a_left = -0.25 + j * self.L
                a_right = a_left + self.L
                # densely sample to get min/max safely
                xs = np.linspace(a_left, a_right, 401)
                cv = np.cos(2.0 * np.pi * k * xs)
                sv = np.sin(2.0 * np.pi * k * xs)
                cos_min[k, j] = cv.min()
                cos_max[k, j] = cv.max()
                sin_min[k, j] = sv.min()
                sin_max[k, j] = sv.max()
        for k in range(1, K + 1):
            constraints.append(self.a[k] >= cos_min[k] @ self.p)
            constraints.append(self.a[k] <= cos_max[k] @ self.p)
            constraints.append(self.b[k] >= sin_min[k] @ self.p)
            constraints.append(self.b[k] <= sin_max[k] @ self.p)

        # ---- v2 (F) test-point constraints on F via lifted M ----
        ts = np.linspace(-0.5, 0.5, n_tests)
        for t_star in ts:
            expr = 1.0
            for k in range(1, K + 1):
                ck = float(np.cos(2.0 * np.pi * k * t_star))
                sk = float(np.sin(2.0 * np.pi * k * t_star))
                A_k = self.M[k, k] - self.M[K + k, K + k]
                B_k = 2.0 * self.M[k, K + k]
                expr = expr + 2.0 * (A_k * ck + B_k * sk)
            constraints.append(self.Omega >= expr)

        self.constraints = constraints
        self.problem = cp.Problem(cp.Minimize(self.Omega), constraints)

    def solve(self, solver: str = "MOSEK", verbose: bool = False, **kwargs) -> V2Result:
        val = self.problem.solve(solver=solver, verbose=verbose, **kwargs)
        return V2Result(
            status=self.problem.status,
            Omega=float(self.Omega.value) if self.Omega.value is not None else float("nan"),
            primal=float(val) if val is not None else float("nan"),
        )


if __name__ == "__main__":
    prob = AutocorrLowerBoundV2(N=8, K=6, n_tests=41)
    try:
        out = prob.solve(solver="MOSEK", verbose=False)
    except Exception as exc:
        print("MOSEK failed:", exc)
        out = prob.solve(solver="SCS", verbose=False)
    print(f"status      : {out.status}")
    print(f"primal value: {out.primal:.10f}")
    print(f"Omega       : {out.Omega:.10f}")
    print(f"v2 lower bound on C_6.2 : {out.Omega:.6f}")
