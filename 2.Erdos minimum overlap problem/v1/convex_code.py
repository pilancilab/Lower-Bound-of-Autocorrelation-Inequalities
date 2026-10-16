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
        # v1 NEW CONSTRAINT (rigorousproof.md, Theorem 1):
        #
        #     sum_{m=1..T} (c_m^2 + d_m^2)^2  <=  2 Omega - 1/2
        #
        # Derivation: 0 <= M <= Omega gives M^2 <= Omega*M, integrating
        # gives int M^2 <= Omega. Parseval on [-2,2] gives
        # int M^2 = 1/4 + 2 sum_k (A_k^2 + B_k^2). For k = 2m, Lemma 3
        # of White yields A_{2m} = -(c_m^2 + d_m^2)/2 (since sin(m*pi)=0
        # and a_{2m}=c_m/2, b_{2m}=d_m/2). Dropping non-negative odd-k
        # and B_k^2 terms gives the inequality above.
        #
        # We use slack variables s_m >= c_m^2 + d_m^2 (epigraph) so that
        # sum s_m^2 <= 2 Omega - 1/2 is DCP-convex (SOC).
        # ============================================================
        self.s_norm = cp.Variable(self.T + 1, nonneg=True, name="s_norm")
        self.constraints += [self.s_norm[0] == 0]
        for m in range(1, self.T + 1):
            self.constraints += [
                cp.square(self.c[m]) + cp.square(self.d[m]) <= self.s_norm[m]
            ]
        self.constraints += [
            cp.sum_squares(self.s_norm[1:self.T + 1]) <= 2.0 * self.Omega - 0.5
        ]

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