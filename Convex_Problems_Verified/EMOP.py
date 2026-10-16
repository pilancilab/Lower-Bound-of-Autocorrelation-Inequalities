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

        # Vectorized a/b linkage to c/d/eps/dele (replaces the 4R-constraint loop).
        # Even modes m = 2,4,...,2R: a[m] approx 0.5*c[m/2], b[m] approx 0.5*d[m/2]
        # Odd modes  m = 1,3,...,2R-1: a[m] approx eps[r] + K_a(m)*(const(m) + ctoa_row @ c)
        #                              b[m] approx dele[r] + (4/pi)*(dtob_row @ d)
        R = self.R
        even_m = np.arange(2, 2 * R + 1, 2)
        even_r = even_m // 2
        odd_m  = np.arange(1, 2 * R + 1, 2)
        odd_r  = (odd_m + 1) // 2

        # Even modes
        self.constraints += [
            cp.abs(self.a[even_m] - 0.5 * self.c[even_r]) <= EPS_EQUALITY,
            cp.abs(self.b[even_m] - 0.5 * self.d[even_r]) <= EPS_EQUALITY,
        ]

        # Odd modes
        K_a_odd     = (2.0 * odd_m / self.pi) * np.sin(self.pi * odd_m / 2.0)  # (R,)
        const_a_odd = 0.5 / (odd_m * odd_m)                                     # (R,)
        ctoa_odd    = self.ctoa[odd_m, 1:self.T + 1]                            # (R, T)
        dtob_odd    = self.dtob[odd_m, 1:self.T + 1]                            # (R, T)
        self.constraints += [
            cp.abs(
                self.a[odd_m] - self.eps[odd_r]
                - (K_a_odd * const_a_odd)
                - cp.multiply(K_a_odd, ctoa_odd @ self.c[1:self.T + 1])
            ) <= EPS_EQUALITY,
            cp.abs(
                self.b[odd_m] - self.dele[odd_r]
                - (4.0 / self.pi) * (dtob_odd @ self.d[1:self.T + 1])
            ) <= EPS_EQUALITY,
        ]

        # Vectorized: Omega upper bounds (replaces 2N-constraint loop)
        self.constraints += [
            self.Omega <= 1,
            self.w[1:self.N + 1] <= self.Omega,
            self.v[1:self.N + 1] <= self.Omega,
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

        # Vectorized alm/bep/bem constraints (replaces the 6R-constraint loop).
        # For each m in 1..2R:
        #   FcosLB(m) = (L/2) * sum_j alm[j,m] * (w_j+v_j)
        #             <= (4*sin(pi*m/2)/(pi*m))*a[m] - 2*(a[m]^2 + b[m]^2)
        #   FsinLB(m) <= -(4*sin/(pi*m))*b[m]
        #   FsinUB(m) >= -(4*sin/(pi*m))*b[m]
        # NB: paper page 12 writes coefficient 8, but that's a typo: (3.6) gives
        # B_m = -(4/(m*pi))*sin(m*pi/2)*b_m, and Lemma 5 bounds the LHS by B_m.
        # Correct coefficient is 4, not 8.
        M_total = 2 * R
        m_full = np.arange(1, M_total + 1)
        sin_full = np.sin(self.pi * m_full / 2.0)
        K_cos = 4.0 * sin_full / (self.pi * m_full)                            # (2R,)

        alm_mat = self.alm[1:self.N + 1, 1:M_total + 1].T                      # (2R, N)
        bep_mat = self.bep[1:self.N + 1, 1:M_total + 1].T                      # (2R, N)
        bem_mat = self.bem[1:self.N + 1, 1:M_total + 1].T                      # (2R, N)

        FcosLB_vec = 0.5 * self.L * (alm_mat @ wv)                             # (2R,)
        FsinLB_vec = 0.5 * self.L * (bem_mat @ w - bep_mat @ v)                # (2R,)
        FsinUB_vec = 0.5 * self.L * (bep_mat @ w - bem_mat @ v)                # (2R,)

        a_slice = self.a[1:M_total + 1]
        b_slice = self.b[1:M_total + 1]
        K_sin = -K_cos
        self.constraints += [
            FcosLB_vec
            <= cp.multiply(K_cos, a_slice)
            - 2.0 * (cp.square(a_slice) + cp.square(b_slice)),
            FsinLB_vec <= cp.multiply(K_sin, b_slice),
            FsinUB_vec >= cp.multiply(K_sin, b_slice),
        ]

        # Vectorized eps/dele bounds (replaces the 4R-constraint loop).
        m_R = np.arange(1, R + 1)
        mm = 2.0 * m_R - 1.0
        denom = 4.0 - (mm / self.T) ** 2
        eps_bounds = (1.0 / self.pi) * (1.0 / denom) * 2.0 * mm * (6.0 * self.T**3) ** (-0.5)
        dele_bounds = (4.0 / self.pi) * (1.0 / denom) * (2.0 * self.T) ** (-0.5)
        self.constraints += [
            cp.abs(self.eps[1:R + 1])  <= eps_bounds,
            cp.abs(self.dele[1:R + 1]) <= dele_bounds,
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

        # Bochner-Herglotz Toeplitz PSD block (Theorem 6 / rigorousproof.md):
        # T_f >= 0  and  I - T_f >= 0,  where T_f is built from (a_k - i*b_k)/2.
        R1 = self.R + 1
        I_S = np.eye(2 * R1)
        ReTf_rows = []
        ImTf_rows = []
        for k in range(R1):
            re_row = []
            im_row = []
            for l in range(R1):
                if k == l:
                    re_row.append(0.25)
                    im_row.append(0.0)
                else:
                    diff = abs(k - l)
                    re_row.append(0.5 * self.a[diff])
                    if k > l:
                        im_row.append(-0.5 * self.b[diff])
                    else:
                        im_row.append(0.5 * self.b[diff])
            ReTf_rows.append(re_row)
            ImTf_rows.append(im_row)
        ReTf = cp.bmat(ReTf_rows)
        ImTf = cp.bmat(ImTf_rows)
        Sf_top = cp.hstack([ReTf, -ImTf])
        Sf_bot = cp.hstack([ImTf, ReTf])
        Sf = cp.vstack([Sf_top, Sf_bot])
        self.constraints += [Sf >> 0]
        self.constraints += [I_S - Sf >> 0]

        self.problem = cp.Problem(cp.Minimize(self.Omega), self.constraints)


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