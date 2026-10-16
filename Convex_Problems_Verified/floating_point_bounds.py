from mpmath import iv, mpf, mp
 
def get_ubd_primal_emop(N, L, T, R, h_1, h_2, p_1, p_2, q_1, q_2,
                        eps_equality=1e-8, dps=20):
    """Returns a certified float X >= sup ||x||_1, rounded up."""
    iv.dps = dps; mp.dps = dps
    Ni, Ti, Ri = iv.mpf(N), iv.mpf(T), iv.mpf(R)
    EPS = iv.mpf(eps_equality)
 
    # Omega block
    X = iv.mpf(1)
 
    # w,v block: from mass constraint L*sum(w+v) <= 1+EPS, L=2/N
    #   sum(w+v) <= (1+EPS)/L = (1+EPS)*N/2
    # (user-specified form N/2 * (1+EPS) ; this is the certified version)
    X = X + (Ni/2) * (1 + EPS)
 
    # c,d block: |c_k|,|d_k| <= 2/pi, T nonzero entries each
    X = X + (4*Ti) / iv.pi
 
    # a,b block: user uses 8R (conservative; tight bound is 2R)
    X = X + 8*Ri
 
    # eps,del block: explicit bounds, summed with interval arithmetic
    for m in range(1, R+1):
        mm = iv.mpf(2*m - 1)
        denom = 4 - (mm/Ti)**2                      # > 0 for valid params
        e_m = (2*mm) / (iv.pi * denom * iv.sqrt(6*Ti**3))
        d_m = iv.mpf(4) / (iv.pi * denom * iv.sqrt(2*Ti))
        X = X + e_m + d_m
 
    # auxiliary (canonicalization) blocks:
    #   square epigraphs (4R of them, u >= a_m^2 or b_m^2): each bounded by 4
    #   sum_squares epigraphs (2 of them, u >= sum c_k^2 or sum d_k^2):
    #       each <= T*(2/pi)^2 = 4T/pi^2  (from |c_k|,|d_k| <= 2/pi over T terms)
    #   abs slacks (4R+1 of them): each <= EPS
    X = X + 4*Ri * 4                                # 4R square epigraphs, bound 4 each
    X = X + 2 * (4*Ti) / iv.pi**2                   # 2 sum_squares epigraphs, 4T/pi^2 each
    X = X + (4*Ri + 1) * EPS                        # abs-slack auxiliaries
 
    # return the UPPER endpoint -> certified upper bound (outward rounded)
    return float(X.b)

def get_ubd_primal_autocorr(N, K, Omega_ub=2, dps=50):
    iv.dps = dps; mp.dps = dps
    Ki = iv.mpf(K)
    Oub = iv.mpf(Omega_ub)
    omega = Oub - 1
 
    X = iv.mpf(0)
 
    # Omega
    X = X + Oub
    # p (simplex)
    X = X + iv.mpf(1)
    # a : K+2  (a_0=1 plus K+1 entries each <=1)
    X = X + (Ki + 2)
    # b : K+1  (b_0=0)
    X = X + (Ki + 1)
    # svec(M): (K+1) tr(M), tr(M) <= 1 + sqrt(K(K+1)omega/2)
    trM = 1 + iv.sqrt(Ki * (Ki + 1) * omega / 2)
    X = X + (Ki + 1) * trM
    # v : (K+1)omega/2
    v_l1 = (Ki + 1) * omega / 2
    X = X + v_l1
    # epigraph aux : <= ||v||_1
    X = X + v_l1
    # herm(Q): (1 + sqrt2*K/2) omega
    X = X + (1 + iv.sqrt(2) * Ki / 2) * omega
 
    return float(X.b)
 
 
# ----------------------------------------------------------------------
# reshape SCS PSD svec block -> symmetric matrix of Fractions
#   SCS: column-major lower triangle; off-diagonals scaled by sqrt2.
#   To recover Y: diagonal entries as-is, off-diagonals divided by sqrt2.
#   In EXACT arithmetic we cannot divide by sqrt2 (irrational). Instead we test
#   PD-ness of the SCALED matrix Ytil = D^{1/2} ... -- but simpler & still exact:
#   Y >= 0  <=>  S Y S >= 0 for any invertible diagonal-like S. The svec scaling is
#   Ytil = T(Y) where off-diags are sqrt2*Y_ij. Ytil is NOT similar to Y, so we
#   cannot avoid the sqrt2. SOLUTION: scale the WHOLE test by 2. Define
#   Z with Z_ii = 2*Y_ii (=2*v_ii) and Z_ij = sqrt2*2/sqrt2... -> messy.
#   Cleanest exact route: multiply each off-diagonal of the reshaped (sqrt2-removed)
#   matrix is irrational, so we instead verify PD of  Yhat := congruence with
#   diag scaling by sqrt2 on rows/cols is not rational either.
#   => Use the symmetric matrix  W where we KEEP the svec scaling but apply the
#   exact congruence:  Y >= 0  <=>  Ytil' >= 0 where Ytil' has Ytil'_ij =
#   (svec val) for offdiag and (svec val) for diag, then test PD of  Ytil' under
#   the inner product with weights -- equivalently test PD of the matrix
#   B with B_ii = Y_ii, B_ij = Y_ij. Since dividing offdiag by sqrt2 is a
#   CONGRUENCE-FREE scaling that changes eigenvalues, we must do it.
#   PRACTICAL EXACT FIX: represent sqrt2 offdiag removal by working with the
#   matrix 2*Y: (2Y)_ii = 2*Y_ii (rational), (2Y)_ij = 2*Y_ij = sqrt2*(svecval)
#   still irrational. Hence pure-rational reshape is impossible from sqrt2-svec.
#
#   Therefore: we reshape using a HIGH-PRECISION interval for 1/sqrt2 and run an
#   INTERVAL LDL^T (rigorous), OR we note Y>=0 iff the matrix with offdiag =
#   svecval/sqrt2 is PSD and test that with interval Cholesky. We use interval
#   LDL^T below as the rigorous PSD test (still no false positives).
# ----------------------------------------------------------------------
 
def _reshape_psd_float(blk, n):
    """Reshape SCS column-major lower-tri svec block (with sqrt2 offdiag scaling)
    into an n x n symmetric float64 matrix Y (off-diagonals divided by sqrt2)."""
    import numpy as np
    inv_sqrt2 = 1.0 / np.sqrt(2.0)
    Y = np.zeros((n, n), dtype=np.float64)
    idx = 0
    for j in range(n):
        for i in range(j, n):
            val = float(blk[idx]); idx += 1
            if i == j:
                Y[i, j] = val
            else:
                Y[i, j] = val * inv_sqrt2
                Y[j, i] = val * inv_sqrt2
    return Y
 
 
def _certified_chol_pd(Y, sigma=0.0, dps=80):
    """
    Rigorous SPD certificate for a symmetric float64 matrix Y.

    Proves:
        lambda_min(Y) >= certified_margin.

    Uses:
        Y - sigma I = LL^T + E

    together with

        lambda_min(Y)
        >= sigma + lambda_min(LL^T) - ||E||_2.

    The bound
        lambda_min(LL^T) >= 1 / ||L^{-1}||_2^2
    is certified via interval Gershgorin on
        C = L^{-T} L^{-1}.

    Returns:
        (is_spd, certified_margin)
    """
    import numpy as np
    from mpmath import iv, mp

    iv.dps = dps
    mp.dps = dps

    Y = np.asarray(Y, dtype=np.float64)
    n = Y.shape[0]

    A = Y - sigma * np.eye(n)

    # Numerical Cholesky
    try:
        Lf = np.linalg.cholesky(A)
    except np.linalg.LinAlgError:
        return False, None

    # Interval enclosure of L
    L = [[iv.mpf(float(Lf[i, j])) for j in range(n)] for i in range(n)]

    # Interval inverse B = L^{-1}
    B = [[iv.mpf(0) for _ in range(n)] for _ in range(n)]

    for i in range(n):
        B[i][i] = 1 / L[i][i]
        for j in range(i):
            s = iv.mpf(0)
            for k in range(j, i):
                s += L[i][k] * B[k][j]
            B[i][j] = -s / L[i][i]

    def abs_ub(x):
        return max(abs(x.a), abs(x.b))

    # C = B^T B
    C = [[iv.mpf(0) for _ in range(n)] for _ in range(n)]

    for i in range(n):
        for j in range(n):
            s = iv.mpf(0)
            for k in range(n):
                s += B[k][i] * B[k][j]
            C[i][j] = s

    # Gershgorin upper bound on lambda_max(C)
    lambda_max_C_ub = 0.0

    for i in range(n):
        row_ub = float(C[i][i].b)
        for j in range(n):
            if i != j:
                row_ub += abs_ub(C[i][j])
        lambda_max_C_ub = max(lambda_max_C_ub, row_ub)

    # Therefore:
    # lambda_min(LL^T) >= 1 / lambda_max(C)
    lam_LL_lb = iv.mpf(1.0 / lambda_max_C_ub)

    # Residual E = A - LL^T
    E_fro2 = iv.mpf(0)

    for i in range(n):
        for j in range(n):

            aij = iv.mpf(float(A[i, j]))
            llij = iv.mpf(0)
            for k in range(n):
                llij += L[i][k] * L[j][k]
            eij = aij - llij
            E_fro2 += abs_ub(eij) ** 2

    E_fro_ub = iv.sqrt(E_fro2)
    # Weyl inequality
    margin = iv.mpf(float(sigma)) + lam_LL_lb - E_fro_ub
    return margin.a > 0, float(margin.a)
 
 
# ----------------------------------------------------------------------
# SOC rigorous check via interval arithmetic
# ----------------------------------------------------------------------
def _soc_ok(block, dps=60):
    """block = (t, z_1, ..., z_{k-1}). Rigorously check t >= ||z||_2.
    Returns (ok, margin) with margin = t - upper_bound(||z||) (rigorous lower bound
    on the slack)."""
    iv.dps = dps; mp.dps = dps
    t = iv.mpf(float(block[0]))
    s = iv.mpf(0)
    for zi in block[1:]:
        zz = iv.mpf(float(zi))
        s = s + zz * zz
    nrm = iv.sqrt(s)                    # interval enclosing ||z||
    slack = t - nrm                     # interval enclosing t - ||z||
    return (slack.a > 0), float(slack.a)
 
 
# ----------------------------------------------------------------------
# main entry
# ----------------------------------------------------------------------
def floating_point_bounds(y, dims, sigma_psd=0.0, dps=60):
    """Verify y in K* rigorously, block by block. dims has .zero, .nonneg, .soc
    (list of SOC sizes), .psd (list of PSD matrix dimensions).
 
    sigma_psd: the PSD margin used when solving the dual (Method 1). The PSD test
    factors Y - sigma_psd I and certifies Y > 0 iff sigma_psd exceeds the rounding
    budget c_n u ||Y||. If sigma_psd=0 (raw witness) the PSD blocks will generally
    fail by the rounding budget -- expected; pass the Method-1 margin to certify."""
    import numpy as np
    y = np.asarray(y).ravel()
    results = []
    overall_ok = True
    off = 0
 
    # zero cone: free, skip
    if dims.zero:
        results.append(("zero", 0, True, float("inf")))
        off += dims.zero
 
    # nonneg: y_i >= 0  (exact float comparison; the float IS the witness)
    if dims.nonneg:
        blk = y[off:off + dims.nonneg]
        mn = float(blk.min())
        ok = mn >= 0.0
        results.append(("nonneg", 0, ok, mn))
        overall_ok = overall_ok and ok
        off += dims.nonneg
 
    # SOC blocks: t >= upper_bound(||z||) via interval arithmetic
    for qi, q in enumerate(dims.soc):
        blk = y[off:off + q]
        ok, margin = _soc_ok(blk, dps=dps)
        results.append(("soc", qi, ok, margin))
        overall_ok = overall_ok and ok
        off += q
 
    # PSD blocks: certified floating-point Cholesky
    for si, n in enumerate(dims.psd):
        ln = n * (n + 1) // 2
        blk = y[off:off + ln]
        Y = _reshape_psd_float(blk, n)
        ok, cert_margin = _certified_chol_pd(Y, sigma_psd, dps=dps)
        results.append(("psd", si, ok, cert_margin))
        overall_ok = overall_ok and ok
        off += ln
 
    return {"ok": overall_ok, "blocks": results, "consumed": off, "total": len(y)}

"""
Measure the floating-point rounding error eps_A in the canonical coefficients
of the two programs, by computing every coefficient in (a) double precision
exactly as the solver code does, and (b) high precision via mpmath, then taking
the maximum absolute difference.
 
The high-precision value is the ground truth: at dps digits, mpmath's own
rounding is ~10^{-dps}, far below double's ~1e-16, so (double - highprec) is the
double rounding error to all relevant precision.
"""
from __future__ import annotations
import math
from math import ceil, floor
import numpy as np
import mpmath as mp
 
mp.mp.dps = 50  # 50 decimal digits of working precision (ground truth)
U = 2.0 ** -53   # IEEE double unit roundoff
 
# ======================================================================
#  C_6.5  (White overlap program)  --  coefficient families
# ======================================================================
def overlap_coeff_errors(N, T, R):
    """Return dict: family -> (max_abs_err, max_err_over_u)."""
    Ld = 2.0 / N
    Lh = mp.mpf(2) / N
    errs = {}
    def rec(fam, e):
        e = mp.mpf(e)
        if e > errs.get(fam, mp.mpf(0)):
            errs[fam] = e
 
    # ---- alpha/beta envelopes  (mirrors _build_alpha_beta) ----
    # double: center = pi*m*L*(j-0.5)/2 ; rad = pi*m*L/4 ; cos/sin(center) +/- rad
    for m in range(1, 2 * R + 1):
        for j in range(1, N + 1):
            c_d = math.pi * m * Ld * (j - 0.5) / 2.0
            r_d = math.pi * m * Ld / 4.0
            c_h = mp.pi * m * Lh * (mp.mpf(j) - mp.mpf(1) / 2) / 2
            r_h = mp.pi * m * Lh / 4
            rec("alp", abs(mp.mpf(math.cos(c_d) + r_d) - (mp.cos(c_h) + r_h)))
            rec("alm", abs(mp.mpf(math.cos(c_d) - r_d) - (mp.cos(c_h) - r_h)))
            rec("bep", abs(mp.mpf(math.sin(c_d) + r_d) - (mp.sin(c_h) + r_h)))
            rec("bem", abs(mp.mpf(math.sin(c_d) - r_d) - (mp.sin(c_h) - r_h)))
 
    # ---- ctoa  (mirrors _build_ctoa):  (-1)^k / (m^2 - 4 k^2) ----
    for m in range(1, 2 * R + 1, 2):  # odd modes only (as used in code)
        for k in range(1, T + 1):
            d_d = ((-1.0) ** k) / (m * m - 4.0 * k * k)
            d_h = (mp.mpf(-1) ** k) / (mp.mpf(m) ** 2 - 4 * mp.mpf(k) ** 2)
            rec("ctoa", abs(mp.mpf(d_d) - d_h))
 
    # ---- dtob  (mirrors _build_dtob): k*(-1)^k*sin(pi m/2)/(m^2-4k^2) ----
    for m in range(1, 2 * R + 1, 2):  # odd modes only
        for k in range(1, T + 1):
            d_d = k * ((-1.0) ** k) * math.sin(math.pi * m / 2.0) / (m * m - 4.0 * k * k)
            d_h = k * (mp.mpf(-1) ** k) * mp.sin(mp.pi * m / 2) / (mp.mpf(m) ** 2 - 4 * mp.mpf(k) ** 2)
            rec("dtob", abs(mp.mpf(d_d) - d_h))
 
    # ---- moment-row coefficients: L^2 j, L^2 (j-1), L^3 j^2, L^3 (j-1)^2 ----
    for j in range(1, N + 1):
        rec("mom", abs(mp.mpf(Ld ** 2 * j) - Lh ** 2 * j))
        rec("mom", abs(mp.mpf(Ld ** 2 * (j - 1)) - Lh ** 2 * (j - 1)))
        rec("mom", abs(mp.mpf(Ld ** 3 * j * j) - Lh ** 3 * j * j))
        rec("mom", abs(mp.mpf(Ld ** 3 * (j - 1) * (j - 1)) - Lh ** 3 * (j - 1) * (j - 1)))
 
    # ---- K_cos = 4 sin(pi m/2)/(pi m) ; mass row L ; eps/dele bounds ----
    for m in range(1, 2 * R + 1):
        kd = 4.0 * math.sin(math.pi * m / 2.0) / (math.pi * m)
        kh = 4 * mp.sin(mp.pi * m / 2) / (mp.pi * m)
        rec("Kcos", abs(mp.mpf(kd) - kh))
    rec("mass", abs(mp.mpf(Ld) - Lh))
    for m in range(1, R + 1):
        mm = 2.0 * m - 1.0
        mmh = 2 * mp.mpf(m) - 1
        den_d = 4.0 - (mm / T) ** 2
        den_h = 4 - (mmh / T) ** 2
        e_d = (1.0 / math.pi) * (1.0 / den_d) * 2.0 * mm * (6.0 * T ** 3) ** -0.5
        e_h = (1 / mp.pi) * (1 / den_h) * 2 * mmh * (6 * mp.mpf(T) ** 3) ** mp.mpf(-0.5)
        rec("eps_bnd", abs(mp.mpf(e_d) - e_h))
        dl_d = (4.0 / math.pi) * (1.0 / den_d) * (2.0 * T) ** -0.5
        dl_h = (4 / mp.pi) * (1 / den_h) * (2 * mp.mpf(T)) ** mp.mpf(-0.5)
        rec("del_bnd", abs(mp.mpf(dl_d) - dl_h))
    return errs
 
 
# ======================================================================
#  C_6.2  (autocorrelation program) -- analytical cell extrema + nu/Q
# ======================================================================
def cos_extrema_double(k, a_l, a_r):
    cand = [math.cos(2*math.pi*k*a_l), math.cos(2*math.pi*k*a_r)]
    m_lo = int(ceil(2.0*k*a_l + 1e-12)); m_hi = int(floor(2.0*k*a_r - 1e-12))
    for m in range(m_lo, m_hi+1):
        x = m/(2.0*k)
        if a_l < x < a_r: cand.append((-1.0)**m)
    return min(cand), max(cand)
 
def cos_extrema_high(k, a_l, a_r):
    cand = [mp.cos(2*mp.pi*k*a_l), mp.cos(2*mp.pi*k*a_r)]
    m_lo = int(mp.ceil(2*k*a_l + mp.mpf('1e-12'))); m_hi = int(mp.floor(2*k*a_r - mp.mpf('1e-12')))
    for m in range(m_lo, m_hi+1):
        x = mp.mpf(m)/(2*k)
        if a_l < x < a_r: cand.append(mp.mpf(-1)**m)
    return min(cand), max(cand)
 
def sin_extrema_double(k, a_l, a_r):
    cand = [math.sin(2*math.pi*k*a_l), math.sin(2*math.pi*k*a_r)]
    m_lo = int(ceil(2.0*k*a_l - 0.5 + 1e-12)); m_hi = int(floor(2.0*k*a_r - 0.5 - 1e-12))
    for m in range(m_lo, m_hi+1):
        x = (2*m+1)/(4.0*k)
        if a_l < x < a_r: cand.append((-1.0)**m)
    return min(cand), max(cand)
 
def sin_extrema_high(k, a_l, a_r):
    cand = [mp.sin(2*mp.pi*k*a_l), mp.sin(2*mp.pi*k*a_r)]
    m_lo = int(mp.ceil(2*k*a_l - mp.mpf('0.5') + mp.mpf('1e-12'))); m_hi = int(mp.floor(2*k*a_r - mp.mpf('0.5') - mp.mpf('1e-12')))
    for m in range(m_lo, m_hi+1):
        x = (2*m+1)/(mp.mpf(4)*k)
        if a_l < x < a_r: cand.append(mp.mpf(-1)**m)
    return min(cand), max(cand)
 
def autocorr_coeff_errors(N, K):
    Kext = K + 1
    dim_p = 2 * N
    Lf_d = 1.0 / (4 * N)
    Lf_h = mp.mpf(1) / (4 * N)
    errs = {}
    def rec(fam, e):
        e = mp.mpf(e)
        if e > errs.get(fam, mp.mpf(0)): errs[fam] = e
 
    # cell extrema for the envelope (cos_min/max, sin_min/max)
    for k in range(1, Kext + 1):
        for j in range(dim_p):
            al_d = -0.25 + j * Lf_d;  ar_d = al_d + Lf_d
            al_h = mp.mpf(-1)/4 + j * Lf_h; ar_h = al_h + Lf_h
            cmn_d, cmx_d = cos_extrema_double(k, al_d, ar_d)
            cmn_h, cmx_h = cos_extrema_high(k, al_h, ar_h)
            smn_d, smx_d = sin_extrema_double(k, al_d, ar_d)
            smn_h, smx_h = sin_extrema_high(k, al_h, ar_h)
            rec("cos_min", abs(mp.mpf(cmn_d) - cmn_h)); rec("cos_max", abs(mp.mpf(cmx_d) - cmx_h))
            rec("sin_min", abs(mp.mpf(smn_d) - smn_h)); rec("sin_max", abs(mp.mpf(smx_d) - smx_h))
 
    # Fejer weights 1 - k/(K+1)
    for k in range(1, K + 1):
        wd = 1.0 - k / (K + 1)
        wh = 1 - mp.mpf(k) / (K + 1)
        rec("weights", abs(mp.mpf(wd) - wh))
    return errs
 
 
def report(title, errs):
    print(f"\n===== {title} =====")
    overall = mp.mpf(0)
    for fam, e in sorted(errs.items(), key=lambda kv: -kv[1]):
        print(f"  {fam:10s}: {float(e):.4e}  = {float(e/U):8.3f} * u")
        overall = max(overall, e)
    print(f"  {'-'*40}")
    print(f"  eps_A (max over all families) = {float(overall):.4e} = {float(overall/U):.3f} * u")
    return overall