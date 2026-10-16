"""
certified_wrapper.py
====================

Final wrapper: go from high-level problem parameters to a *certified* lower
bound on the optimal value, by

    high-level params
        -> build the conic (CVXPY) model
        -> canonicalize to SCS standard form  (A, b, c, dims)
        -> compute the matching a-priori primal bound  X >= ||x*||_1
        -> verify the SCS dual witness  y  and assemble the rigorous bound L
           (delegated to certified_lower_bound, Method 5)

It is deliberately agnostic about *how* the conic model is produced.  You hand
it a builder callable (params -> cp.Problem) for each problem class, and it does
the rest.  Two thin convenience entry points are provided:

    certified_lower_bound_emop(...)        # EMOP:  N,L,T,R,h,p,q
    certified_lower_bound_autocorr(...)    # Autocorr 6.2:  N,K

Both call the same engine, `run_certified`.

ASSUMPTIONS (correct these if your setup differs):
  * The conic problems are modeled in CVXPY and solved with SCS.
  * `certified_lower_bound`, `get_ubd_primal_emop`, `get_ubd_primal_autocorr`,
    and `floating_point_bounds` are importable (your existing modules).
  * The objective is a *minimization* whose optimal value is Omega*; the dual
    gives the certified *lower* bound.  (Method 5's sign convention: it computes
    -b^T y with the SCS-internal b, which is what CVXPY hands back.)

If instead you already have raw (A, b, c, dims) on disk, skip the builder and
call `run_certified(..., problem=None, scs_data=(A,b,c,dims), ...)`.
"""

from __future__ import annotations

import numpy as np
import cvxpy as cp

# --- your existing modules -------------------------------------------------
from certified_bound import certified_lower_bound           # Method 5 (doc 1)
from floating_point_bounds import (                          # X bounds  (doc 2)
    get_ubd_primal_emop,
    get_ubd_primal_autocorr,
)


# ---------------------------------------------------------------------------
# 1.  Extract SCS standard-form data from a CVXPY problem
# ---------------------------------------------------------------------------
def scs_data_from_problem(problem: cp.Problem):
    """Canonicalize a CVXPY problem to SCS standard form.

    Returns (A, b, c, dims) where dims is the SCS ConeDims object exposing
    .zero, .nonneg, .soc, .psd in the order Method 2/5 expect.

    CONVENTION NOTE (this is the part that silently breaks certificates):
    CVXPY's SCS data describes the primal as
            minimize  c^T x   s.t.   A x + s = b,   s in K.
    The dual is        maximize -b^T y  s.t.  A^T y + c = 0,   y in K*,
    which is EXACTLY Method 5's form (residual A^T y + c, objective -b^T y).
    So the (A, b, c) from get_problem_data, used verbatim, is the matched
    system for the dual witness.  We do NOT transpose, negate, or reorder; the
    whole point is that y must belong to THIS A,b,c.  The dims cone order
    (zero, nonneg, soc, psd) also matches what Methods 1/2 walk.

    Method 1/2 only handle zero/nonneg/soc/psd, so we assert no exp/pow blocks.
    """
    data, _, _ = problem.get_problem_data(cp.SCS)
    A = data["A"]
    b = data["b"]
    c = data["c"]
    dims = data["dims"]

    # guard: Method 1/2 don't model exponential / power cones
    exp = getattr(dims, "exp", 0)
    p3 = getattr(dims, "p", []) or []
    if exp or len(p3) > 0:
        raise NotImplementedError(
            f"SCS data contains exp({exp}) / power({len(p3)}) cones, which the "
            "certified pipeline (Methods 1,2) does not handle."
        )

    # shape sanity: A is (m, n); b in R^m; c in R^n; cone sizes sum to m
    m, n = A.shape
    if b.shape[0] != m or c.shape[0] != n:
        raise ValueError(
            f"Shape mismatch from SCS data: A={A.shape}, b={b.shape}, c={c.shape}"
        )
    cone_total = (dims.zero + dims.nonneg
                  + sum(dims.soc) + sum(s * (s + 1) // 2 for s in dims.psd))
    if cone_total != m:
        raise ValueError(
            f"Cone dims sum to {cone_total} but A has {m} rows; the certified "
            "block walk in Methods 1/2 would be misaligned."
        )
    return A, b, c, dims


# ---------------------------------------------------------------------------
# 2.  On the dual witness y
# ---------------------------------------------------------------------------
# We deliberately do NOT reuse SCS's returned dual. Two reasons:
#   (a) SCS's optimal dual sits on the boundary of K*, so the certified PSD
#       Cholesky test (_certified_chol_pd) fails by exactly the rounding budget
#       -- it needs an interior point with a margin.
#   (b) Reusing it would require reconstructing CVXPY's internal scaling/reorder
#       to map the reported dual back onto the extracted (A,b,c); any mismatch
#       silently certifies the wrong system.
# So we pass y=None to certified_lower_bound, which invokes Method 1
# (solve_dual_with_margin) on the SAME (A,b,c,dims) we extracted -> matched by
# construction. The residual ||A^T y + c||_inf computed inside Method 5 is the
# witness's own consistency check; we surface it in `info` and assert on it.


# ---------------------------------------------------------------------------
# 3.  The engine
# ---------------------------------------------------------------------------
def run_certified(
    *,
    problem: cp.Problem | None,
    X: float,
    scs_data=None,
    eps_A=0.0, eps_b=0.0, eps_c=0.0,
    margin=1e-6, sigma_psd=None, dps=60,
    solver_opts=None,
    residual_tol=1e-6,
    y_save_path=None,
    solver_eps=1e-12,
):
    """Core driver.  Either pass a CVXPY `problem` (it will be solved and
    canonicalized) or pass `scs_data=(A,b,c,dims)` directly.

    Parameters
    ----------
    problem    : cp.Problem to solve & canonicalize, or None if scs_data given.
    X          : certified primal bound X >= ||x*||_1 (from get_ubd_primal_*).
    scs_data   : optional (A,b,c,dims) tuple bypassing the CVXPY route.
    eps_A/b/c  : rigorous UPPER bounds on data error; 0 certifies the float
                 problem exactly.
    margin     : dual interior margin handed to Method 1.
    sigma_psd  : PSD factoring shift (defaults inside certified_lower_bound).
    dps        : mpmath working precision.
    solver_opts: dict of extra kwargs for problem.solve() (e.g. eps, max_iters).
    residual_tol : if the dual-feasibility residual ||A^T y + c||_inf exceeds
                 this, the witness does not satisfy the equality constraints
                 well enough for the bound to be meaningful -> we flag it.
                 (It does not by itself invalidate the rigorous L, which already
                 folds the residual into the penalty, but a large residual means
                 a loose / useless bound and usually signals a solver problem.)

    Returns (L, valid, info).  `valid` is the cone certificate from Method 2;
    `info["residual_ok"]` separately reports the feasibility-residual check.
    """
    if scs_data is not None:
        A, b, c, dims = scs_data
    else:
        if problem is None:
            raise ValueError("Provide either `problem` or `scs_data`.")
        # Skip the float-problem solve entirely: scs_data_from_problem only needs
        # CVXPY canonicalization (get_problem_data forces it), and we don't
        # reuse SCS's returned primal/dual anyway. Saves the primal SCS solve.
        A, b, c, dims = scs_data_from_problem(problem)

    # y=None -> Method 1 solves an interior witness on THIS (A,b,c,dims).
    # If y_save_path is given, the witness is saved as .npz for later re-audit.
    L, valid, info = certified_lower_bound(
        A, b, c, dims, X,
        eps_A=eps_A, eps_b=eps_b, eps_c=eps_c,
        margin=margin, sigma_psd=sigma_psd, dps=dps, y=None,
        y_save_path=y_save_path, solver_eps=solver_eps,
    )

    info = dict(info)
    # the residual is computed inside Method 5 and reported as residual_float_ub
    res = info.get("residual_float_ub", float("inf"))
    info["residual_ok"] = bool(res <= residual_tol)
    info["residual_tol"] = float(residual_tol)
    info["X_used"] = float(X)
    info["dims"] = {
        "zero": dims.zero, "nonneg": dims.nonneg,
        "soc": list(dims.soc), "psd": list(dims.psd),
    }
    info["A_shape"] = tuple(A.shape)
    return L, valid, info


# ---------------------------------------------------------------------------
# 4.  Problem-specific entry points
# ---------------------------------------------------------------------------
def certified_lower_bound_emop(
    N, L, T, R, h_1, h_2, p_1, p_2, q_1, q_2,
    *,
    build_emop,                       # callable: params -> cp.Problem
    eps_equality=1e-8,
    eps_A=1e-12, eps_b=1e-12, eps_c=1e-12,
    margin=1e-8, sigma_psd=None,
    dps_X=20, dps_cert=60,
    solver_opts=None, residual_tol=1e-6,
    y_save_path=None,
):
    """EMOP: build the conic model from (N,L,T,R,h,p,q), compute the matching
    X via get_ubd_primal_emop, and certify.

    `build_emop` is YOUR formulation function with signature
        build_emop(N, L, T, R, h_1, h_2, p_1, p_2, q_1, q_2) -> cp.Problem
    (the function that creates exactly the w,v / c,d / a,b / eps,del blocks and
    the canonicalization epigraphs that get_ubd_primal_emop accounts for).
    """
    X = get_ubd_primal_emop(
        N, L, T, R, h_1, h_2, p_1, p_2, q_1, q_2,
        eps_equality=eps_equality, dps=dps_X,
    )
    problem = build_emop(N, L, T, R, h_1, h_2, p_1, p_2, q_1, q_2)
    return run_certified(
        problem=problem, X=X,
        eps_A=eps_A, eps_b=eps_b, eps_c=eps_c,
        margin=margin, sigma_psd=sigma_psd, dps=dps_cert,
        solver_opts=solver_opts, residual_tol=residual_tol,
        y_save_path=y_save_path,
    )


def certified_lower_bound_autocorr(
    N, K,
    *,
    build_autocorr,                   # callable: params -> cp.Problem
    Omega_ub=2,
    eps_A=1e-12, eps_b=1e-12, eps_c=1e-12,
    margin=1e-8, sigma_psd=None,
    dps_X=50, dps_cert=60,
    solver_opts=None, residual_tol=1e-6,
    y_save_path=None,
    solver_eps=1e-12,
):
    """Autocorrelation 6.2: build the moment/SDP relaxation from (N,K), compute
    X via get_ubd_primal_autocorr, and certify.

    `build_autocorr` is YOUR formulation function with signature
        build_autocorr(N, K) -> cp.Problem
    producing the Omega / p / a / b / svec(M) / v / herm(Q) blocks that
    get_ubd_primal_autocorr accounts for.
    """
    X = get_ubd_primal_autocorr(N, K, Omega_ub=Omega_ub, dps=dps_X)
    problem = build_autocorr(N, K)
    return run_certified(
        problem=problem, X=X,
        eps_A=eps_A, eps_b=eps_b, eps_c=eps_c,
        margin=margin, sigma_psd=sigma_psd, dps=dps_cert,
        solver_opts=solver_opts, residual_tol=residual_tol,
        y_save_path=y_save_path, solver_eps=solver_eps,
    )
