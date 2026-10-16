import numpy as np
import cvxpy as cp
from mpmath import iv, mp
 
from floating_point_bounds import floating_point_bounds   # Method 2
 
 
# --- rigorous directed-rounded helpers -------------------------------------
def _dot_round_down(u, v):
    """LOWER bound on true u.v (u,v float, taken exact)."""
    s = iv.mpf(0)
    for ui, vi in zip(u, v):
        s = s + iv.mpf(float(ui)) * iv.mpf(float(vi))
    return s.a
 
 
def _l1_round_up(v):
    """UPPER bound on ||v||_1."""
    s = iv.mpf(0)
    for vi in v:
        s = s + abs(iv.mpf(float(vi)))
    return s.b
 
 
def _residual_inf_round_up(A, y, c):
    """UPPER bound on ||A^T y + c||_inf (float-problem residual), per column."""
    Acsc = A.tocsc() if hasattr(A, "tocsc") else None
    n = A.shape[1]
    worst = iv.mpf(0)
    if Acsc is not None:
        indptr, indices, dat = Acsc.indptr, Acsc.indices, Acsc.data
        for j in range(n):
            acc = iv.mpf(float(c[j]))
            for p in range(indptr[j], indptr[j + 1]):
                acc = acc + iv.mpf(float(dat[p])) * iv.mpf(float(y[indices[p]]))
            mag = max(abs(acc.a), abs(acc.b))
            if mag > worst:
                worst = mag
    else:
        Ad = np.asarray(A)
        for j in range(n):
            acc = iv.mpf(float(c[j]))
            col = Ad[:, j]
            for i in np.nonzero(col)[0]:
                acc = acc + iv.mpf(float(col[i])) * iv.mpf(float(y[i]))
            mag = max(abs(acc.a), abs(acc.b))
            if mag > worst:
                worst = mag
    return worst
 
 
# --- Method 1 : interior witness via margin dual solve ---------------------
def solve_dual_with_margin(A, b, c, dims, margin=1e-6, solver_eps=1e-12,
                            mosek_threads=8, y_save_path=None):
    """max -b^T y  s.t.  A^T y + c = 0,  y interior to K* by `margin`.
    Cone order (SCS): zero(free), nonneg, SOC, PSD.

    Uses MOSEK interior-point (multi-threaded, 5-20x faster than SCS for these
    SDPs). `solver_eps` controls MOSEK's PFEAS/DFEAS tolerance; REL_GAP/INFEAS
    are set 10x tighter so the cone audit at margin>=solver_eps passes cleanly."""
    import math
    m = A.shape[0]
    y = cp.Variable(m)
    cons = []
    idx = dims.zero
    if dims.nonneg > 0:
        cons.append(y[idx:idx + dims.nonneg] >= margin)
        idx += dims.nonneg
    for q in dims.soc:
        yi = y[idx:idx + q]
        cons.append(cp.SOC(yi[0] - margin, yi[1:]))
        idx += q
    for s in dims.psd:
        sz = s * (s + 1) // 2
        yi = y[idx:idx + sz]
        rows = [[None] * s for _ in range(s)]
        k = 0
        for col in range(s):
            for r in range(col, s):
                if r == col:
                    rows[r][col] = yi[k]
                else:
                    rows[r][col] = yi[k] / math.sqrt(2)
                    rows[col][r] = yi[k] / math.sqrt(2)
                k += 1
        cons.append(cp.bmat(rows) >> margin * np.eye(s))
        idx += sz
    cons.append(A.T @ y + c == 0)
    mosek_params = {
        "MSK_IPAR_NUM_THREADS": int(mosek_threads),
        "MSK_DPAR_INTPNT_CO_TOL_DFEAS":  solver_eps,
        "MSK_DPAR_INTPNT_CO_TOL_PFEAS":  solver_eps,
        "MSK_DPAR_INTPNT_CO_TOL_REL_GAP": solver_eps * 0.1,
        "MSK_DPAR_INTPNT_CO_TOL_INFEAS": solver_eps * 0.1,
        "MSK_DPAR_INTPNT_CO_TOL_MU_RED": solver_eps * 0.1,
    }
    prob = cp.Problem(cp.Maximize(-(b @ y)), cons)
    prob.solve(solver=cp.MOSEK, mosek_params=mosek_params)
    # Detect cell-excluded (primal infeasible): the original primal
    #   min c^T x  s.t.  Ax + s = b,  s in K
    # has no admissible point, so its optimum is +inf and the dual is
    # unbounded. CVXPY/MOSEK report this as status='unbounded'.
    if prob.status in ("unbounded", "unbounded_inaccurate"):
        return None  # sentinel: cell is vacuously verified (Omega = +inf)
    if y.value is None:
        raise RuntimeError(
            f"MOSEK dual solve failed: status={prob.status!r}, "
            f"value={prob.value!r}")
    y_arr = np.asarray(y.value).ravel()
    if y_save_path is not None:
        import os
        os.makedirs(os.path.dirname(y_save_path) or ".", exist_ok=True)
        np.savez_compressed(y_save_path if y_save_path.endswith(".npz") else y_save_path + ".npz",
                            y=y_arr, margin=margin, solver_eps=solver_eps,
                            dims_zero=dims.zero, dims_nonneg=dims.nonneg,
                            dims_soc=np.array(list(dims.soc), dtype=np.int64),
                            dims_psd=np.array(list(dims.psd), dtype=np.int64))
    return y_arr
 
 
# --- Method 5 : assemble certified L (uses Methods 1, 2, 3) ----------------
def certified_lower_bound(A, b, c, dims, X,
                          eps_A=0.0, eps_b=0.0, eps_c=0.0,
                          margin=1e-6, sigma_psd=None, dps=60, y=None,
                          y_save_path=None, solver_eps=1e-12):
    """Rigorous lower bound L on Omega*. Returns (L, valid, info).
 
    A,b,c,dims : SCS standard-form float data + cone dims.
    X          : a-priori bound X >= ||x*||_1 (Get_Ubd_primal).
    eps_A,eps_b,eps_c : rigorous UPPER bounds on data error (Get_eABC/proof);
                 pass 0 to certify the float problem only.
    margin     : dual interior margin (Method 1).
    sigma_psd  : PSD factoring shift (defaults to margin).
    y          : optional precomputed witness; if None Method 1 solves it.
    """
    iv.dps = dps; mp.dps = dps
    if sigma_psd is None:
        # factor Y - sigma I with sigma strictly INSIDE the witness interior
        # (Method 1 made lambda_min(Y) >~ margin; leave headroom so Cholesky of
        #  Y - sigma I succeeds while sigma still beats the rounding budget).
        sigma_psd = margin / 2.0
 
    # Method 1: interior witness
    if y is None:
        y = solve_dual_with_margin(A, b, c, dims, margin=margin,
                                    solver_eps=solver_eps,
                                    y_save_path=y_save_path)

    # Cell-excluded short-circuit: solve_dual_with_margin returns None when
    # the dual is unbounded, which happens iff the primal min c^T x s.t. ...
    # is +infinity (no admissible x). Vacuously verified: Omega = +inf.
    if y is None:
        info = {
            "cone": None, "residual_float_ub": 0.0,
            "y_l1_ub": 0.0, "eps_total_ub": 0.0,
            "obj_down": float("inf"), "penalty_ub": 0.0,
            "X": float(X), "L": float("inf"), "valid": True,
            "method": "primal_infeasible_cell_excluded",
        }
        return float("inf"), True, info

    # Method 2: verify y in K*
    cone = floating_point_bounds(y, dims, sigma_psd=sigma_psd, dps=dps)
 
    # Method 3: residual + data-error, all upper-bounded
    res_float = _residual_inf_round_up(A, y, c)
    y_l1 = _l1_round_up(y)
    eps_total = (iv.mpf(float(res_float))
                 + iv.mpf(float(eps_A)) * iv.mpf(float(y_l1))
                 + iv.mpf(float(eps_c)))
 
    # objective rounded DOWN, penalty rounded UP, assemble
    obj_down = _dot_round_down([-bi for bi in b], y)
    penalty = iv.mpf(eps_total.b) * iv.mpf(float(X)) + iv.mpf(float(eps_b)) * iv.mpf(float(y_l1))
    L = float((iv.mpf(obj_down) - penalty).a)
 
    info = {"cone": cone, "residual_float_ub": float(res_float),
            "y_l1_ub": float(y_l1), "eps_total_ub": float(eps_total.b),
            "obj_down": float(obj_down), "penalty_ub": float(penalty.b),
            "X": float(X), "L": L, "valid": cone["ok"]}
    return L, cone["ok"], info
 