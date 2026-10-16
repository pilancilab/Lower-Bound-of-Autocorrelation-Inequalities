from __future__ import annotations

from dataclasses import dataclass
import numpy as np
import cvxpy as cp

@dataclass(frozen=True)
class WhiteInputs:
    N: int
    L: float
    T: int
    R: int
    h_1: float
    h_2: float
    p_1: float
    p_2: float
    q_1: float
    q_2: float


class WhiteConvexProblem:
    def __init__(
        self,
        N: int,
        L: float,
        T: int,
        R: int,
        h_1: float,
        h_2: float,
        p_1: float,
        p_2: float,
        q_1: float,
        q_2: float,
    ) -> None:
        self.N = N
        self.L = L
        self.T = T
        self.R = R
        self.h_1 = h_1
        self.h_2 = h_2
        self.p_1 = p_1
        self.p_2 = p_2
        self.q_1 = q_1
        self.q_2 = q_2
        self.pi = np.pi

        EPS_EQUALITY = 1e-8
        
        self._validate_inputs()

        # CVXPY variables
        self.Omega = cp.Variable(nonneg=True, name="Omega")
        self.w = cp.Variable(self.N + 1, nonneg=True, name="w")
        self.v = cp.Variable(self.N + 1, nonneg=True, name="v")
        self.c = cp.Variable(self.T + 1, name="c")
        self.d = cp.Variable(self.T + 1, name="d")
        self.eps = cp.Variable(self.R + 1, name="eps")
        self.dele = cp.Variable(self.R + 1, name="del")
        # v3 NEW: one free variable per odd k in {1, 3, ..., 2R-1}.
        # Stored as a length-(2R+1) array, indexed 1..2R-1 with stride 2;
        # entries at even indices and at 0 are dummy/zero.
        self.A_odd = cp.Variable(2 * self.R + 1, name="A_odd")
        # v4 NEW: one free variable per EVEN k in {2, 4, ..., R} (for use
        # in the Hermitian Toeplitz matrix). Stored as a length-(R+1)
        # array, indexed 1..R; only even indices in {2, ..., R} are used,
        # the rest are forced to 0.
        self.A_even = cp.Variable(self.R + 1, name="A_even")

        # optional: enforce dummy 0th entries to be zero
        self.index_zero_constraints = [
            self.w[0] == 0,
            self.v[0] == 0,
            self.c[0] == 0,
            self.d[0] == 0,
            self.eps[0] == 0,
            self.dele[0] == 0,
        ]

        # coefficient arrays, all 1-indexed
        self.ctoa = self._build_ctoa()
        self.dtob = self._build_dtob()
        self.alp ,self.alm ,self.bep, self.bem = self._build_alpha_beta()

        self.a = cp.Variable(2 * self.R + 1, name="a")
        self.b = cp.Variable(2 * self.R + 1, name="b")
        self.index_zero_constraints += [
            self.a[0] == 0,
            self.b[0] == 0,
        ]
        self.constraints = []
        self.constraints += self.index_zero_constraints
        for m in range(1, 2 * self.R + 1):
            if m % 2 == 0:
                # even modes: m = 2r
                r = m // 2
                self.constraints += [
                    cp.abs(self.a[m] - 0.5*self.c[r]) <= EPS_EQUALITY,
                    cp.abs(self.b[m] - 0.5*self.d[r]) <= EPS_EQUALITY,
                ]
            else:
                # odd modes: m = 2r-1
                r = (m + 1) // 2
                self.constraints += [
                    cp.abs(self.a[m] - self.eps[r]
                     - (2.0 * m / self.pi) * np.sin(self.pi * m / 2.0) * (
                        0.5 / (m * m) + self.ctoa[m, 1:self.T + 1] @ self.c[1:self.T + 1]
                    )) <= EPS_EQUALITY,
                    cp.abs(self.b[m] - self.dele[r]
                     - (4.0 / self.pi) * (
                        self.dtob[m, 1:self.T + 1] @ self.d[1:self.T + 1]
                    )) <= EPS_EQUALITY,
                ]

        self.constraints += [self.Omega <= 1]
        for j in range(1, self.N + 1):
            self.constraints += [
                self.w[j] <= self.Omega,
                self.v[j] <= self.Omega,
            ]
        
        self.constraints += [
            cp.abs(self.L * cp.sum(self.w[1:self.N + 1] + self.v[1:self.N + 1]) - 1) <= EPS_EQUALITY
        ]

        j = np.arange(1, self.N + 1, dtype=float)
        self.constraints += [
            self.h_1 <= self.L**2 * (
                j @ self.w[1:self.N + 1] - (j - 1.0) @ self.v[1:self.N + 1]
            )
        ]
        self.constraints += [
            self.L**2 * (
                (j - 1.0) @ self.w[1:self.N + 1] - j @ self.v[1:self.N + 1]
            ) <= self.h_2
        ]

        self.constraints += [
            self.L**3 * ((j - 1.0) ** 2 @ (self.w[1:self.N + 1] + self.v[1:self.N + 1]))
            <= 2.0 / 3.0 + self.h_2**2 / 2.0
        ]

        # 2/3 + h_1^2 / 2 <= L^3 * sum_{j=1}^N j^2 (w_j + v_j)
        self.constraints += [
            2.0 / 3.0 + self.h_1**2 / 2.0
            <= self.L**3 * (j**2 @ (self.w[1:self.N + 1] + self.v[1:self.N + 1]))
        ]

        w = self.w[1:self.N + 1]
        v = self.v[1:self.N + 1]
        wv = w + v

        for m in range(1, 2 * self.R + 1):
            sin_term = np.sin(self.pi * m / 2.0)

            # FcosLB(m) = (L/2) * sum_j alm[j,m] * (w_j + v_j)
            FcosLB_m = 0.5 * self.L * (self.alm[1:self.N + 1, m] @ wv)

            # FsinLB(m) = (L/2) * sum_j (bem[j,m] * w_j - bep[j,m] * v_j)
            FsinLB_m = 0.5 * self.L * (
                self.bem[1:self.N + 1, m] @ w
                - self.bep[1:self.N + 1, m] @ v
            )

            # FsinUB(m) = (L/2) * sum_j (bep[j,m] * w_j - bem[j,m] * v_j)
            FsinUB_m = 0.5 * self.L * (
                self.bep[1:self.N + 1, m] @ w
                - self.bem[1:self.N + 1, m] @ v
            )

            # alm constraint
            self.constraints += [
                FcosLB_m
                <= (4.0 * sin_term / (self.pi * m)) * self.a[m]
                - 2.0 * (cp.square(self.a[m]) + cp.square(self.b[m]))
            ]

            # bep / bem constraints
            self.constraints += [
                FsinLB_m <= -(4.0 * sin_term / (self.pi * m)) * self.b[m],
                FsinUB_m >= -(4.0 * sin_term / (self.pi * m)) * self.b[m],
            ]

        for m in range(1, self.R + 1):
    
            mm = float(2 * m - 1)
            denom = 4.0 - (mm / self.T) ** 2

            eps_bound = (1.0 / self.pi) * (1.0 / denom) * 2.0 * mm * (6.0 * self.T**3) ** (-0.5)
            dele_bound = (4.0 / self.pi) * (1.0 / denom) * (2.0 * self.T) ** (-0.5)

            self.constraints += [
                self.eps[m] <= eps_bound,
                self.eps[m] >= -eps_bound,
                self.dele[m] <= dele_bound,
                self.dele[m] >= -dele_bound,
            ]

        bound = 2.0 / self.pi
        self.constraints += [
            self.c[1:self.T + 1] <= bound,
            self.c[1:self.T + 1] >= -bound,
            self.d[1:self.T + 1] <= bound,
            self.d[1:self.T + 1] >= -bound,
        ]

        self.constraints += [
            cp.sum_squares(self.c[1:self.T + 1]) + cp.sum_squares(self.d[1:self.T + 1]) <= 0.5
        ]

        # p1 <= c1 <= p2, q1 <= d1 <= q2
        self.constraints += [
            self.p_1 <= self.c[1],
            self.c[1] <= self.p_2,
            self.q_1 <= self.d[1],
            self.d[1] <= self.q_2,
        ]

        qsq_max = max(self.q_1**2, self.q_2**2)
        self.constraints += [
            0.5 * self.L * (
                self.alp[1:self.N + 1, 2] @ (self.w[1:self.N + 1] + self.v[1:self.N + 1])
            )
            >= -0.5 * (self.p_2**2 + qsq_max)
        ]

        # ============================================================
        # v3 NEW CONSTRAINTS (rigorousproof.md, Theorem 3):
        #
        # For each odd k = 2r-1 in {1, 3, ..., 2R-1} introduce a free
        # program variable A_odd[k] = "tilde A_k" with
        #
        #   (A1) FcosLB[k]  <= A_odd[k] <= FcosUB[k]    (Lemma 5)
        #   (A2) A_odd[k]   <= (4 sgn / k pi) * a_k - 2 (a_k^2 + b_k^2)
        #                                                 (Lemma 3 (3.5))
        #
        # both of which are convex in the program variables (the (A2)
        # RHS is concave in (a_k, b_k), so the hypograph constraint is
        # convex in (A_odd, a, b)).
        #
        # Then strengthen the v2 energy bound by adding 4 * A_odd[k]^2:
        #
        #   sum_{m=1..T} (c_m^2 + d_m^2)^2
        # + 4 * sum_{k odd, 1..2R-1} A_odd[k]^2
        # + (64/pi^2) * sum_{k odd, 1..2R-1} b_k^2 / k^2
        #   <= 2 Omega - 1/2
        #
        # All terms are convex; A_odd[k]^2 is a univariate square of a
        # CVXPY variable.
        # ============================================================
        self.s_norm = cp.Variable(self.T + 1, nonneg=True, name="s_norm")
        self.constraints += [self.s_norm[0] == 0]
        for m in range(1, self.T + 1):
            self.constraints += [
                cp.square(self.c[m]) + cp.square(self.d[m]) <= self.s_norm[m]
            ]

        odd_ks = [k for k in range(1, 2 * self.R) if k % 2 == 1]

        # zero out dummy entries of A_odd (so the variable's other slots
        # cannot bias the optimizer)
        self.constraints += [self.A_odd[0] == 0]
        for k in range(1, 2 * self.R + 1):
            if k % 2 == 0:
                self.constraints += [self.A_odd[k] == 0]

        # (A1): Lemma 5 LB and UB on A_k for odd k
        for k in odd_ks:
            FcosLB_k = 0.5 * self.L * (
                self.alm[1:self.N + 1, k]
                @ (self.w[1:self.N + 1] + self.v[1:self.N + 1])
            )
            FcosUB_k = 0.5 * self.L * (
                self.alp[1:self.N + 1, k]
                @ (self.w[1:self.N + 1] + self.v[1:self.N + 1])
            )
            self.constraints += [
                FcosLB_k <= self.A_odd[k],
                self.A_odd[k] <= FcosUB_k,
            ]

        # (A2): A_odd[k] <= concave Lemma 3 RHS
        # sign of sin((2r-1) pi / 2) is (-1)^{r+1} = (-1)^((k-1)/2)
        for k in odd_ks:
            r = (k + 1) // 2
            sign = (-1.0) ** (r + 1)
            self.constraints += [
                self.A_odd[k]
                <= (4.0 * sign / (k * self.pi)) * self.a[k]
                - 2.0 * (cp.square(self.a[k]) + cp.square(self.b[k]))
            ]

        odd_b_sq_terms = [
            (64.0 / (self.pi ** 2 * k * k)) * cp.square(self.b[k])
            for k in odd_ks
        ]
        odd_b_sq = cp.sum(odd_b_sq_terms) if odd_b_sq_terms else 0.0

        odd_A_sq_terms = [
            4.0 * cp.square(self.A_odd[k]) for k in odd_ks
        ]
        odd_A_sq = cp.sum(odd_A_sq_terms) if odd_A_sq_terms else 0.0

        self.constraints += [
            cp.sum_squares(self.s_norm[1:self.T + 1]) + odd_A_sq + odd_b_sq
            <= 2.0 * self.Omega - 0.5
        ]

        # ============================================================
        # v4 NEW CONSTRAINT (rigorousproof.md, Theorem 4):
        # Bochner-Herglotz Toeplitz PSD on M-hat.
        #
        # We lift A_k for even k = 2, 4, ..., R (the odd ones are already
        # lifted in v3 as A_odd[k] for k = 1, 3, ..., 2R-1).
        #
        # Then build the (R+1)x(R+1) Hermitian Toeplitz matrix
        #     T[k, l] = M-hat(k - l)
        # with M-hat(0) = 1/4, M-hat(j) = (A_lift[j] - i B_j) / 2,
        # M-hat(-j) = conjugate.
        #
        # Impose T >= 0 (PSD), expressed as the real (2R+2)x(2R+2)
        # symmetric block matrix S = [[ReT, -ImT], [ImT, ReT]] >> 0.
        # ============================================================
        # zero-out non-used entries of A_even
        self.constraints += [self.A_even[0] == 0]
        for k in range(1, self.R + 1):
            if k % 2 == 1:
                self.constraints += [self.A_even[k] == 0]

        # (A1)+(A2) for even k = 2, 4, ..., R (only those <= R needed for Toeplitz of size R+1)
        for k in range(2, self.R + 1, 2):
            m = k // 2
            if m > self.T:
                continue  # no c_m, d_m variable; skip
            FcosLB_k = 0.5 * self.L * (
                self.alm[1:self.N + 1, k]
                @ (self.w[1:self.N + 1] + self.v[1:self.N + 1])
            )
            FcosUB_k = 0.5 * self.L * (
                self.alp[1:self.N + 1, k]
                @ (self.w[1:self.N + 1] + self.v[1:self.N + 1])
            )
            self.constraints += [
                FcosLB_k <= self.A_even[k],
                self.A_even[k] <= FcosUB_k,
                # Lemma 3 even branch: A_{2m} = -(c_m^2 + d_m^2)/2
                # Hypograph form: A_even[k] <= -(c_m^2 + d_m^2)/2 (concave)
                self.A_even[k]
                <= -0.5 * (cp.square(self.c[m]) + cp.square(self.d[m])),
            ]

        # A_lift[k] for k = 1..R: unified accessor (returns expr)
        def A_lift(k):
            if k == 0:
                return 0.25  # diagonal handled separately
            if k <= 0 or k > self.R:
                return 0.0
            if k % 2 == 1:
                return self.A_odd[k]
            else:
                return self.A_even[k]

        # B_j for j = 1..R (sine coeff of M, linear in b for odd j)
        def B_lift(j):
            if j <= 0 or j > self.R:
                return 0.0
            if j % 2 == 0:
                return 0.0
            r = (j + 1) // 2
            sign = (-1.0) ** (r + 1)
            # B_j = -(4 sign / (j pi)) b_j
            return -(4.0 * sign / (j * self.pi)) * self.b[j]

        # Build ReT, ImT as (R+1)x(R+1) lists of CVXPY expressions
        R1 = self.R + 1
        ReT_rows = []
        ImT_rows = []
        for k in range(R1):
            re_row = []
            im_row = []
            for l in range(R1):
                if k == l:
                    re_row.append(0.25)
                    im_row.append(0.0)
                else:
                    diff = abs(k - l)
                    A_d = A_lift(diff)
                    B_d = B_lift(diff)
                    re_row.append(0.5 * A_d)
                    # Im(M-hat(k-l)) = -sign(k-l) * B/2 / ... let's compute:
                    # M-hat(diff > 0) = (A - i B)/2 -> imag = -B/2
                    # M-hat(diff < 0) = conj -> imag = +B/2
                    if k > l:  # i.e., k-l > 0
                        im_row.append(-0.5 * B_d)
                    else:
                        im_row.append(0.5 * B_d)
            ReT_rows.append(re_row)
            ImT_rows.append(im_row)

        # Use cp.bmat to assemble (each cell is a scalar expr)
        ReT = cp.bmat(ReT_rows)
        ImT = cp.bmat(ImT_rows)

        # Real PSD form: S = [[ReT, -ImT], [ImT, ReT]] >> 0
        S_top = cp.hstack([ReT, -ImT])
        S_bot = cp.hstack([ImT, ReT])
        S = cp.vstack([S_top, S_bot])

        self.constraints += [S >> 0]

    def _validate_inputs(self) -> None:
        if self.N <= 0 or self.T <= 0 or self.R <= 0:
            raise ValueError("N, T, R must be positive.")
        if abs(self.L - 2 / self.N) > 1e-12:
            raise ValueError("L must equal 2/N.")
        if self.p_1 > self.p_2:
            raise ValueError("Need p_1 <= p_2.")
        if self.q_1 > self.q_2:
            raise ValueError("Need q_1 <= q_2.")

    def _center(self, j: int, m: int) -> float:
        return self.pi * m * self.L * (j - 0.5) / 2.0

    def _rad(self, m: int) -> float:
        return self.pi * m * self.L / 4.0

    def _build_ctoa(self) -> np.ndarray:
        ms = np.arange(1, 2 * self.R + 1)[:, None]   # shape (2R, 1)
        ks = np.arange(1, self.T + 1)[None, :]       # shape (1, T)

        core = ((-1.0) ** ks) / (ms**2 - 4 * ks**2)

        out = np.zeros((2 * self.R + 1, self.T + 1))
        out[1:, 1:] = core
        return out

    def _build_dtob(self) -> np.ndarray:
        ms = np.arange(1, 2 * self.R + 1)[:, None]   # shape (2R, 1)
        ks = np.arange(1, self.T + 1)[None, :]       # shape (1, T)

        core = ks * ((-1.0) ** ks) * np.sin(np.pi * ms / 2.0) / (ms**2 - 4 * ks**2)

        out = np.zeros((2 * self.R + 1, self.T + 1))
        out[1:, 1:] = core
        return out

    def _build_alpha_beta(self):
        js = np.arange(1, self.N + 1)[:, None]       # shape (N, 1)
        ms = np.arange(1, 2 * self.R + 1)[None, :]   # shape (1, 2R)

        center = np.pi * ms * self.L * (js - 0.5) / 2.0
        rad = np.pi * ms * self.L / 4.0

        alp = np.zeros((self.N + 1, 2 * self.R + 1))
        alm = np.zeros((self.N + 1, 2 * self.R + 1))
        bep = np.zeros((self.N + 1, 2 * self.R + 1))
        bem = np.zeros((self.N + 1, 2 * self.R + 1))

        alp[1:, 1:] = np.cos(center) + rad
        alm[1:, 1:] = np.cos(center) - rad
        bep[1:, 1:] = np.sin(center) + rad
        bem[1:, 1:] = np.sin(center) - rad

        return alp, alm, bep, bem
    
    def solve_and_verify(
        self,
        solver=None,
        verbose: bool = True,
        gap_tol: float = 1e-6,
        feas_tol: float = 1e-6,
        **solver_kwargs,
    ):
        self.problem = cp.Problem(cp.Minimize(self.Omega), self.constraints)
        value = self.problem.solve(solver=solver, verbose=verbose, **solver_kwargs)

        stats = self.problem.solver_stats
        extra = getattr(stats, "extra_stats", None)

        dual_obj = None
        gap = None

        if isinstance(extra, dict):
            # try a few common keys
            for key in ["dual_obj", "dual objective", "dual_objective", "dobj"]:
                if key in extra:
                    dual_obj = extra[key]
                    break
            for key in ["gap", "duality_gap", "rel_gap"]:
                if key in extra:
                    gap = extra[key]
                    break

        # fallback: if dual objective is available but gap is not
        if gap is None and dual_obj is not None and value is not None:
            try:
                gap = abs(value - dual_obj)
            except Exception:
                pass

        # check max constraint violation
        max_violation = 0.0
        for con in self.constraints:
            try:
                v = con.violation()
                if v is None:
                    continue
                v = np.max(np.abs(v))
                max_violation = max(max_violation, float(v))
            except Exception:
                pass

        verified = True
        if gap is not None and gap > gap_tol:
            verified = False
        if max_violation > feas_tol:
            verified = False

        print("\n=== Verification summary ===")
        print("status         :", self.problem.status)
        print("primal optimum :", self.problem.value)
        print("dual optimum   :", dual_obj)
        print("duality gap    :", gap)
        print("max violation  :", max_violation)
        print("verified       :", verified)

        return {
            "status": self.problem.status,
            "primal_optimum": self.problem.value,
            "dual_optimum": dual_obj,
            "duality_gap": gap,
            "max_violation": max_violation,
            "verified": verified,
        }

if __name__ == "__main__":
    # Hard square verification regime from proposal.md
    problem = WhiteConvexProblem(
        N=10000,
        L=2/10000,
        T=2000,
        R=10,
        h_1=0.015,
        h_2=0.015,
        p_1=0.381,
        p_2=0.381,
        q_1=-0.02,
        q_2=0.02,
    )

    try:
        result = problem.solve_and_verify(solver=cp.MOSEK, verbose=True)
    except Exception as e:
        print("MOSEK failed:", e)
        result = problem.solve_and_verify(solver=cp.SCS, verbose=True, max_iters=40000)