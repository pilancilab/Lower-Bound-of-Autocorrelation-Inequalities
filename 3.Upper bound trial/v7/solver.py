"""
v7 — Refines v6 (0.38091069141515083) by:

  1. Upsampling exactly (each bin repeated factor times) to n = 1200/1800/2400
     WITHOUT touching the L-BFGS / LSE step (which previously kicked the
     solution out of basin). Then long mass-transport polish at the higher n
     where there are more degrees of freedom.
  2. Adding **two-pair coordinated mass moves** (transfer δ from {i1,i2}
     simultaneously to {j1,j2}) — needed when the LP descent direction is
     supported on >2 coords and pair-only moves are stuck at a saddle.
  3. Final repeated LSE/AS/MT cycles.
"""

import json, os, time
import numpy as np
from scipy.optimize import minimize, brentq, linprog

EPS = 1e-9


# ---- shared core (identical to v6) ----------------------------------------

def _L(n):
    L = 1
    while L < 2 * n: L *= 2
    return L


def overlap_full(a):
    n = len(a); b = 1.0 - a; L = _L(n)
    A = np.fft.rfft(a, n=L); B = np.fft.rfft(b, n=L)
    full = np.fft.irfft(np.conjugate(A) * B, n=L)
    M = np.empty(2 * n - 1)
    M[n - 1:] = full[:n]; M[:n - 1] = full[L - (n - 1):L]
    return M * (2.0 / n)


def grad_weighted(a, w):
    n = len(a); b = 1.0 - a; L = _L(n)
    w_circ = np.zeros(L)
    w_circ[:n] = w[n - 1:]; w_circ[L - (n - 1):L] = w[:n - 1]
    W = np.fft.rfft(w_circ, n=L); B = np.fft.rfft(b, n=L); A = np.fft.rfft(a, n=L)
    T1 = np.fft.irfft(np.conjugate(W) * B, n=L)[:n]
    T2 = np.fft.irfft(W * A, n=L)[:n]
    return (2.0 / n) * (T1 - T2)


def grad_single_shift(a, d):
    n = len(a); b = 1.0 - a; g = np.zeros(n)
    if d >= 0:
        g[:n - d] += b[d:]; g[d:] -= a[:n - d]
    else:
        kk = -d; g[kk:] += b[:n - kk]; g[:n - kk] -= a[kk:]
    return g * (2.0 / n)


def _bin_grad(a, p):
    n = len(a); b = 1.0 - a
    g = np.zeros(2 * n - 1)
    idx_b = np.arange(n) - p + (n - 1); g[idx_b] += b
    idx_a = p - np.arange(n) + (n - 1); g[idx_a] -= a
    return g * (2.0 / n)


def sigmoid(z):
    out = np.empty_like(z); pos = z >= 0
    out[pos] = 1.0 / (1.0 + np.exp(-z[pos]))
    ez = np.exp(z[~pos]); out[~pos] = ez / (1.0 + ez)
    return out


def project_sum_c(z, target):
    def F(c): return sigmoid(z + c).sum() - target
    lo, hi = -50.0, 50.0
    while F(lo) > 0 and lo > -1e5: lo -= 50.0
    while F(hi) < 0 and hi < 1e5:  hi += 50.0
    return brentq(F, lo, hi, xtol=1e-13, rtol=1e-13)


def make_fg(n, beta):
    target = n / 2.0
    def fg(z):
        c = project_sum_c(z, target); a = sigmoid(z + c)
        vals = overlap_full(a)
        m = vals.max(); e = np.exp(beta * (vals - m)); Z = e.sum()
        smax = m + np.log(Z) / beta
        w = e / Z
        df_da = grad_weighted(a, w)
        ab = a * (1.0 - a); S = ab.sum()
        if S < 1e-30: return float(smax), np.zeros_like(z)
        proj = (df_da * ab).sum() / S
        return float(smax), ab * (df_da - proj)
    return fg


def run_lbfgs(z, n, betas, iters):
    for beta in betas:
        fg = make_fg(n, beta)
        res = minimize(fg, z, jac=True, method='L-BFGS-B',
                       options={'maxiter': iters, 'ftol': 1e-15,
                                'gtol': 1e-11, 'maxcor': 80})
        z = res.x
    return z


def optimize_from_a(a, betas, iters):
    n = len(a)
    z = np.log(np.clip(a, EPS, 1 - EPS)) - np.log(np.clip(1 - a, EPS, 1 - EPS))
    z = z + project_sum_c(z, n / 2.0)
    z = run_lbfgs(z, n, betas, iters)
    a_out = sigmoid(z + project_sum_c(z, n / 2.0))
    return a_out, float(overlap_full(a_out).max())


# ---- mass transport: pair + two-pair coordinated moves ---------------------

def _try_pair(a, i, j, step, vals, n, cur_ub):
    delta = min(step, a[i], 1.0 - a[j])
    if delta <= 0: return None
    gi = _bin_grad(a, int(i))
    gj = _bin_grad(a, int(j))
    lin_k = gj - gi
    d_vals = delta * lin_k
    k1 = int(j) - int(i)
    d_vals[k1 + n - 1] += (2.0 / n) * delta * delta
    d_vals[-k1 + n - 1] += (2.0 / n) * delta * delta
    new_max = float((vals + d_vals).max())
    return new_max, delta


def _try_two_pair(a, i1, i2, j1, j2, step, vals, n, cur_ub):
    """Move δ from i1 and δ from i2 to j1 and j2 respectively."""
    if len({i1, i2, j1, j2}) != 4: return None
    delta = min(step, a[i1], a[i2], 1.0 - a[j1], 1.0 - a[j2])
    if delta <= 0: return None
    # Compute delta_vals = δ * (g_{j1}+g_{j2}-g_{i1}-g_{i2}) + quadratic terms
    g_j1 = _bin_grad(a, int(j1))
    g_j2 = _bin_grad(a, int(j2))
    g_i1 = _bin_grad(a, int(i1))
    g_i2 = _bin_grad(a, int(i2))
    lin_k = g_j1 + g_j2 - g_i1 - g_i2
    d_vals = delta * lin_k
    # Quadratic: from -(2/n) sum_l v_l v_{l+k} where v has support on {i1,i2}=-δ, {j1,j2}=+δ
    # For each pair (p, q) with v_p v_q nonzero: contribution at k = q - p is
    # -(2/n) (v_p * v_q). Then symmetric.
    coords = [(int(i1), -delta), (int(i2), -delta), (int(j1), delta), (int(j2), delta)]
    for p, vp in coords:
        for q, vq in coords:
            if p == q: continue
            k = q - p
            if abs(k) <= n - 1:
                d_vals[k + n - 1] += -(2.0 / n) * vp * vq  # but each pair counted twice (p,q) and (q,p)
    # Note: my loop counts each (unordered) pair twice for distinct p,q, which matches
    # the formula sum_l v_l v_{l+k} where the sum is over ordered (l, l+k).
    new_max = float((vals + d_vals).max())
    return new_max, delta


def mass_transport_polish(a, time_budget=180.0, step_init=1e-3,
                          step_min=2.0 ** -36, verbose=False, log_every=120,
                          use_two_pair=True, two_pair_every=3):
    a = a.copy(); n = len(a); t_start = time.time()
    step = step_init; accept = 0; stall = 0; itr = 0
    rng = np.random.default_rng(424242)
    while time.time() - t_start < time_budget:
        itr += 1
        vals = overlap_full(a); cur_ub = vals.max()
        active = np.where(vals >= cur_ub - 5e-5)[0]
        if len(active) == 0: break
        w = np.zeros(2 * n - 1); ev = np.exp(5000.0 * (vals[active] - cur_ub))
        w[active] = ev / ev.sum()
        g = grad_weighted(a, w)
        elig_d = np.where(a > step)[0]; elig_r = np.where(a < 1.0 - step)[0]
        if len(elig_d) == 0 or len(elig_r) == 0:
            step *= 0.5
            if step < step_min: break
            continue
        K_top = min(35, len(elig_d), len(elig_r))
        donor_top = elig_d[np.argsort(g[elig_d])[::-1][:K_top]]
        recv_top = elig_r[np.argsort(g[elig_r])[:K_top]]
        donor_rand = rng.choice(elig_d, size=min(15, len(elig_d)), replace=False)
        recv_rand = rng.choice(elig_r, size=min(15, len(elig_r)), replace=False)
        donors = np.unique(np.concatenate([donor_top, donor_rand]))
        recvs = np.unique(np.concatenate([recv_top, recv_rand]))

        best_drop = 0.0; best_move = None  # ('pair', i, j, delta) or ('two', i1, i2, j1, j2, delta)
        # Try single pairs
        for i in donors:
            for j in recvs:
                if i == j: continue
                ret = _try_pair(a, int(i), int(j), step, vals, n, cur_ub)
                if ret is None: continue
                new_max, delta = ret
                drop = cur_ub - new_max
                if drop > best_drop + 1e-17:
                    best_drop = drop
                    best_move = ('pair', int(i), int(j), delta)
        # Periodically try two-pair moves restricted to the top few donors/receivers
        if use_two_pair and (itr % two_pair_every == 0) and best_drop < 1e-12:
            top_donors = donor_top[:8] if len(donor_top) > 0 else []
            top_recvs = recv_top[:8] if len(recv_top) > 0 else []
            for ii1, i1 in enumerate(top_donors):
                for i2 in top_donors[ii1 + 1:]:
                    for jj1, j1 in enumerate(top_recvs):
                        for j2 in top_recvs[jj1 + 1:]:
                            ret = _try_two_pair(a, int(i1), int(i2),
                                                int(j1), int(j2),
                                                step, vals, n, cur_ub)
                            if ret is None: continue
                            new_max, delta = ret
                            drop = cur_ub - new_max
                            if drop > best_drop + 1e-17:
                                best_drop = drop
                                best_move = ('two', int(i1), int(i2),
                                             int(j1), int(j2), delta)
        if best_move is None:
            stall += 1
            if stall >= 4:
                step *= 0.5; stall = 0
                if step < step_min: break
            continue
        if best_move[0] == 'pair':
            _, i, j, delta = best_move
            a[i] -= delta; a[i] = max(0.0, a[i])
            a[j] += delta; a[j] = min(1.0, a[j])
        else:
            _, i1, i2, j1, j2, delta = best_move
            a[i1] -= delta; a[i1] = max(0.0, a[i1])
            a[i2] -= delta; a[i2] = max(0.0, a[i2])
            a[j1] += delta; a[j1] = min(1.0, a[j1])
            a[j2] += delta; a[j2] = min(1.0, a[j2])
        accept += 1; stall = 0
        if verbose and accept % log_every == 0:
            print(f"  MT it {itr} acc {accept} step={step:.2e} "
                  f"ub={float(overlap_full(a).max()):.17f} "
                  f"drop={best_drop:.2e} kind={best_move[0]}")
    return a


def active_set_polish(a, max_iters=300, tol=1e-14, active_tol=1e-4,
                      max_active=400):
    a = a.copy(); n = len(a)
    last_ub = float(overlap_full(a).max())
    for it in range(max_iters):
        vals = overlap_full(a); ub = vals.max()
        active_idx = np.where(vals >= ub - active_tol)[0]
        if len(active_idx) > max_active:
            order = np.argsort(vals[active_idx])[::-1][:max_active]
            active_idx = active_idx[order]
        K = len(active_idx)
        if K == 0: break
        grads = np.zeros((K, n))
        for j, idx in enumerate(active_idx):
            grads[j] = grad_single_shift(a, int(idx - (n - 1)))
        c = np.zeros(n + 1); c[-1] = 1.0
        A_ub = np.concatenate([grads, -np.ones((K, 1))], axis=1); b_ub = np.zeros(K)
        A_eq = np.concatenate([np.ones((1, n)), np.zeros((1, 1))], axis=1); b_eq = np.zeros(1)
        bounds_lp = [(-1.0, 1.0)] * n + [(None, None)]
        try:
            res = linprog(c, A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=b_eq,
                          bounds=bounds_lp, method='highs')
        except Exception: break
        if not res.success or res.x is None: break
        d = res.x[:n]; t = res.x[-1]
        if t > -1e-14: break
        with np.errstate(divide='ignore', invalid='ignore'):
            pos = d > 1e-15; neg = d < -1e-15
            step_pos = np.min((1.0 - a)[pos] / d[pos]) if pos.any() else np.inf
            step_neg = np.min(a[neg] / (-d[neg])) if neg.any() else np.inf
        max_step = 0.999 * min(step_pos, step_neg)
        if max_step <= 1e-16: break
        L = _L(n); D = np.fft.rfft(d, n=L)
        ac = np.fft.irfft(np.conjugate(D) * D, n=L)
        q = np.empty(2 * n - 1)
        q[n - 1:] = ac[:n]; q[:n - 1] = ac[L - (n - 1):L]
        q *= -(2.0 / n)
        A_fft = np.fft.rfft(a, n=L); B_fft = np.fft.rfft(1 - a, n=L)
        t1 = np.fft.irfft(np.conjugate(D) * B_fft, n=L)
        t2 = np.fft.irfft(np.conjugate(A_fft) * D, n=L)
        lin = np.empty(2 * n - 1)
        lin[n - 1:] = (t1 - t2)[:n]; lin[:n - 1] = (t1 - t2)[L - (n - 1):L]
        lin *= (2.0 / n)
        s = max_step; best_s, best_u = 0.0, ub
        for _ in range(80):
            u = float((vals + s * lin + s * s * q).max())
            if u < best_u - 1e-15: best_u, best_s = u, s; break
            s *= 0.5
            if s < 1e-17: break
        if best_s > 0:
            for f in (1.3, 1.5, 2.0, 0.8, 0.6):
                ss = min(best_s * f, max_step)
                u = float((vals + ss * lin + ss * ss * q).max())
                if u < best_u - 1e-15: best_u, best_s = u, ss
        if best_s == 0.0: break
        a = a + best_s * d; a = np.clip(a, 0.0, 1.0)
        excess = a.sum() - n / 2.0
        if abs(excess) > 1e-12:
            inner = (a > 1e-6) & (a < 1 - 1e-6)
            if inner.any():
                a[inner] -= excess / inner.sum(); a = np.clip(a, 0.0, 1.0)
        cur_ub = float(overlap_full(a).max())
        if last_ub - cur_ub < tol: last_ub = cur_ub; break
        last_ub = cur_ub
    return a


def upsample_repeat(a, factor):
    return np.repeat(a, factor)


# ---- driver ----------------------------------------------------------------

def solve(verbose=False, time_budget_sec=900):
    here = os.path.dirname(os.path.abspath(__file__))
    t_start = time.time()

    v6_path = os.path.join(here, "..", "v6", "sequence.json")
    with open(v6_path) as f:
        v6_data = json.load(f)
    a_v6 = np.asarray(v6_data["a"], dtype=float)
    ub_v6 = float(v6_data.get("upper_bound", float(overlap_full(a_v6).max())))
    if verbose:
        print(f"[warm v6] ub={ub_v6:.17f}, n={len(a_v6)}")

    best_a, best_ub = a_v6, ub_v6

    # ---- Phase 1: upsample 2x to n=1200, polish with MT (no LSE in between) ----
    if verbose: print(f"\n[Phase 1: n=1200 MT polish]")
    a_1200 = upsample_repeat(a_v6, 2).astype(np.float64)
    ub_up = float(overlap_full(a_1200).max())
    if verbose: print(f"   upsample U = {ub_up:.17f}")
    mt_budget = time_budget_sec * 0.30
    a_mt = mass_transport_polish(a_1200, time_budget=mt_budget,
                                 step_init=2.0 ** -10, step_min=2.0 ** -36,
                                 verbose=verbose, use_two_pair=True)
    ub_mt = float(overlap_full(a_mt).max())
    if verbose: print(f"   n=1200 MT done: {ub_mt:.17f}")
    if ub_mt < best_ub:
        if verbose: print(f"   (NEW BEST)")
        best_a, best_ub = a_mt, ub_mt

    # ---- Phase 2: upsample 3x to n=1800, polish ----
    if time.time() - t_start < time_budget_sec * 0.50:
        if verbose: print(f"\n[Phase 2: n=1800 MT polish]")
        a_1800 = upsample_repeat(a_v6, 3).astype(np.float64)
        ub_up = float(overlap_full(a_1800).max())
        if verbose: print(f"   upsample U = {ub_up:.17f}")
        mt_budget = time_budget_sec * 0.20
        a_mt = mass_transport_polish(a_1800, time_budget=mt_budget,
                                     step_init=2.0 ** -10, step_min=2.0 ** -36,
                                     verbose=verbose, use_two_pair=True)
        ub_mt = float(overlap_full(a_mt).max())
        if verbose: print(f"   n=1800 MT done: {ub_mt:.17f}")
        if ub_mt < best_ub:
            if verbose: print(f"   (NEW BEST)")
            best_a, best_ub = a_mt, ub_mt

    # ---- Phase 3: try the best so far back at n=600 (downsample mean) ----
    n_cur = len(best_a)
    if n_cur > 600 and n_cur % 600 == 0:
        factor = n_cur // 600
        a_600 = best_a.reshape(600, factor).mean(axis=1)
        ub_600 = float(overlap_full(a_600).max())
        if verbose: print(f"\n[Phase 3: downsample to n=600] ub={ub_600:.17f}")
        # MT polish at n=600
        a_polished = mass_transport_polish(a_600, time_budget=60.0,
                                           step_init=1e-5,
                                           step_min=2.0 ** -36,
                                           verbose=verbose,
                                           use_two_pair=True)
        ub_p = float(overlap_full(a_polished).max())
        if ub_p < best_ub:
            if verbose: print(f"   n=600 MT done: {ub_p:.17f} (NEW)")
            best_a, best_ub = a_polished, ub_p

    # ---- Phase 4: final MT polish at chosen n ----
    n_cur = len(best_a)
    remaining = max(60.0, time_budget_sec - (time.time() - t_start))
    if verbose:
        print(f"\n[Phase 4: final MT at n={n_cur}] budget {remaining:.0f}s, "
              f"current {best_ub:.17f}")
    a_mt = mass_transport_polish(best_a, time_budget=remaining,
                                 step_init=2.0 ** -16, step_min=2.0 ** -38,
                                 verbose=verbose, use_two_pair=True,
                                 two_pair_every=2)
    ub_mt = float(overlap_full(a_mt).max())
    if ub_mt < best_ub:
        if verbose: print(f"   final MT: {ub_mt:.17f} (NEW)")
        best_a, best_ub = a_mt, ub_mt

    n = len(best_a)
    vals = overlap_full(best_a)
    arg_d = int(np.argmax(vals) - (n - 1))
    seq_path = os.path.join(here, "sequence.json")
    with open(seq_path, "w") as f:
        json.dump({"n": n, "a": best_a.tolist(),
                   "upper_bound": float(best_ub), "argmax_shift": arg_d}, f)
    return {
        "upper_bound": float(best_ub),
        "status": "ok",
        "sequence_path": seq_path,
        "n": int(n),
        "mass_error": float(abs(best_a.sum() - n / 2.0)),
        "box_violation": float(max(0.0, best_a.max() - 1.0, -best_a.min())),
        "argmax_shift": arg_d,
    }


if __name__ == "__main__":
    budget = int(os.environ.get("BUDGET_SEC", "900"))
    print(json.dumps(solve(verbose=True, time_budget_sec=budget), indent=2))
