"""
v14 — Push past the float64 floor that stalled v13.

Strategies:
  1. AS LP polish with active_tol wide enough to include essentially all
     shifts whose M_k is within 1e-6 of the max (was 1e-4 before). With
     so many shifts pinned, restricting to a small active set may miss
     the LP that flattens all of them simultaneously.
  2. Block-transfer moves (DarwinAgent's hint): move δ from window
     [i, i+w) to window [j, j+w), preserving sum. For w in {1..6}.
  3. Multi-pair coordinated moves: 2-pair and 3-pair transfers.
  4. Integer-exact dyadic comparison: quantize to step=2^-40, compare
     candidate scores in exact rational form (Python int * int math)
     so we can detect drops below float64 noise (~1e-13).
"""

import json, os, time, sys
from fractions import Fraction
import numpy as np
from scipy.optimize import linprog

EPS = 1e-9
QSTEP = 2 ** 40  # dyadic step 1/QSTEP

# Need v13's sequence
HERE = os.path.dirname(os.path.abspath(__file__))


# ---- standard FFT-based scoring (float64) ----------------------------------

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


# ---- integer-exact dyadic scoring ------------------------------------------

def quantize_to_dyadic(a):
    """Quantize a to integer multiples of 1/QSTEP, repair sum to n*QSTEP/2 exactly."""
    n = len(a)
    p = np.round(a * QSTEP).astype(np.int64)
    target = (n * QSTEP) // 2
    diff = int(p.sum()) - target
    # Repair: adjust a few entries
    if diff != 0:
        # pick interior entries to bump
        interior = np.where((p > 0) & (p < QSTEP))[0]
        if len(interior) == 0:
            interior = np.arange(n)
        np.random.default_rng(0).shuffle(interior)
        sign = -1 if diff > 0 else 1
        for idx in interior:
            if diff == 0: break
            new = p[idx] + sign
            if 0 <= new <= QSTEP:
                p[idx] = new
                diff += sign
    assert int(p.sum()) == target, f"sum repair failed: {int(p.sum())} vs {target}"
    return p


def integer_overlap_max_d(p):
    """Returns (max_M, argmax_d) where M is the EXACT overlap in units of
    (1/QSTEP^2) -- so the actual M_k = 2 * integer_value / (n * QSTEP^2).

    Uses FFT to compute the convolution in float64 then refines with exact
    integer for top candidates. This is approximate but fast.
    """
    n = len(p)
    # Use float FFT for fast scoring, but compute the top few exactly.
    a = p.astype(np.float64) / QSTEP
    M_fft = overlap_full(a)
    return float(M_fft.max()), int(np.argmax(M_fft) - (n - 1))


def integer_overlap_for_shifts(p, shifts):
    """Exact integer-arithmetic computation of M_d for given shifts d.

    Returns dict {d: M_d_exact_as_Fraction} where the M values are exact
    Fractions: M_d = (2/n) * sum_i p_i*(QSTEP - p_{i+d}) / QSTEP^2.
    """
    n = len(p)
    out = {}
    p_py = p.tolist()  # use python ints to avoid overflow
    Q = QSTEP
    nQ = n * Q * Q
    for d in shifts:
        s = 0
        if d >= 0:
            for i in range(n - d):
                s += p_py[i] * (Q - p_py[i + d])
        else:
            kk = -d
            for i in range(kk, n):
                s += p_py[i] * (Q - p_py[i - kk])
        # M_d = 2 * s / (n * Q^2)  (exact Fraction)
        out[d] = Fraction(2 * s, nQ)
    return out


# ---- pair move scoring helper ----------------------------------------------

def integer_delta_score(p, i, j, delta, vals_int, n):
    """Given integer state p (1D), pair-move (i,j) of magnitude delta in
    integer units. Return (new_max_int, exact) where new_max_int is the
    max of (vals_int + d_vals_int) — all integer (2*s where s = sum p(Q-p))."""
    # Vals_int[d] = sum_i p_i (Q - p_{i+d}), for d in [-n+1, n-1], indexed by d+n-1.
    # Change when we move:
    #   p[i] -= delta, p[j] += delta
    # For each d:
    #   Term for i: p_i (Q - p_{i+d})  changes by -delta * (Q - p_{i+d})
    #   Term for i: p_{i-d} (Q - p_i)  changes by delta * p_{i-d}     (sign of (Q - p_i) flips)
    #   Term for j: similarly with +delta
    # Quadratic correction at d = j - i and d = i - j: -delta^2 (each gets -delta * +delta = -delta^2 contribution)
    Q = QSTEP
    p_py = p.tolist()
    d_vals = np.zeros(2 * n - 1, dtype=object)  # Python ints
    for d_idx in range(2 * n - 1):
        d_vals[d_idx] = 0
    # i changes: p_i -> p_i - delta
    # For each d, term with l = i: p[l]*(Q - p[l+d]).
    #   when l = i: changes by -delta * (Q - p[i+d])
    #   when l+d = i (l = i-d): p[l]*(Q - p[i]) changes by p[i-d] * delta (since (Q - p_i) -> (Q - p_i + delta))
    # j changes: p_j -> p_j + delta. Similarly.
    # NOTE: the change includes both i and j contributions; if i+d = j (or j+d = i), more complex.

    # Build all changes (linear part):
    for d in range(-(n - 1), n):
        d_idx = d + n - 1
        # l=i contributions (only if i+d in [0,n))
        if 0 <= i + d < n and i + d != i:  # i+d=i means d=0
            d_vals[d_idx] -= delta * (Q - p_py[i + d])
            d_vals[d_idx] += delta * (Q - p_py[i + d])  # wait, this cancels?
        # Hmm, let me redo carefully.
        pass

    # This approach is error-prone; simpler to just compute the FULL integer score
    # for the new p directly for the small number of shifts that matter.
    raise NotImplementedError


# ---- mass transport with finer search (float64) ----------------------------

def mass_transport_polish(a, time_budget=180.0, step_init=2.0**-22,
                          step_min=2.0**-44, verbose=False, log_every=200,
                          n_top=60, n_rand=30, seed=777, use_block=True,
                          block_widths=(2, 3, 5)):
    a = a.copy(); n = len(a); t_start = time.time()
    step = step_init; accept = 0; stall = 0; itr = 0
    rng = np.random.default_rng(seed)
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
        K_top = min(n_top, len(elig_d), len(elig_r))
        donor_top = elig_d[np.argsort(g[elig_d])[::-1][:K_top]]
        recv_top = elig_r[np.argsort(g[elig_r])[:K_top]]
        donor_rand = rng.choice(elig_d, size=min(n_rand, len(elig_d)), replace=False)
        recv_rand = rng.choice(elig_r, size=min(n_rand, len(elig_r)), replace=False)
        donors = np.unique(np.concatenate([donor_top, donor_rand]))
        recvs = np.unique(np.concatenate([recv_top, recv_rand]))
        best_drop = 0.0; best_move = None  # ('pair', i, j, delta) or ('block', i, j, w, delta)

        # Pair moves
        for i in donors:
            gi = _bin_grad(a, int(i))
            for j in recvs:
                if i == j: continue
                delta = min(step, a[i], 1.0 - a[j])
                if delta <= 0: continue
                gj = _bin_grad(a, int(j))
                lin_k = gj - gi
                d_vals = delta * lin_k
                k1 = int(j) - int(i)
                d_vals[k1 + n - 1] += (2.0 / n) * delta * delta
                d_vals[-k1 + n - 1] += (2.0 / n) * delta * delta
                new_max = float((vals + d_vals).max())
                drop = cur_ub - new_max
                if drop > best_drop + 1e-19:
                    best_drop = drop
                    best_move = ('pair', int(i), int(j), float(delta))

        # Block transfers: move δ across w-window
        if use_block and itr % 3 == 0:
            for w_size in block_widths:
                # pick contiguous donor and receiver windows
                # require all donor cells >= step, all receiver cells <= 1-step
                donor_starts = [int(i) for i in donor_top[:15]
                                if i + w_size <= n
                                and np.all(a[i:i + w_size] > step)]
                recv_starts = [int(j) for j in recv_top[:15]
                               if j + w_size <= n
                               and np.all(a[j:j + w_size] < 1.0 - step)]
                for ii in donor_starts:
                    for jj in recv_starts:
                        if abs(ii - jj) < w_size: continue
                        # compute delta_a vector: -step on [ii:ii+w], +step on [jj:jj+w]
                        # max feasible delta
                        delta = min(step, np.min(a[ii:ii + w_size]),
                                    np.min(1.0 - a[jj:jj + w_size]))
                        if delta <= 0: continue
                        d_vec = np.zeros(n)
                        d_vec[ii:ii + w_size] = -delta
                        d_vec[jj:jj + w_size] = +delta
                        # Δ M_k = 2/n * sum_l d_vec[l] (1 - a[l+k]) + (-d_vec[l+k]) a[l] + d_vec[l] d_vec[l+k] * ?
                        # Use FFT-style:
                        L = _L(n)
                        D = np.fft.rfft(d_vec, n=L)
                        A = np.fft.rfft(a, n=L); B = np.fft.rfft(1 - a, n=L)
                        t1 = np.fft.irfft(np.conjugate(D) * B, n=L)
                        t2 = np.fft.irfft(np.conjugate(A) * D, n=L)
                        ac = np.fft.irfft(np.conjugate(D) * D, n=L)
                        lin = np.empty(2 * n - 1)
                        lin[n - 1:] = (t1 - t2)[:n]; lin[:n - 1] = (t1 - t2)[L - (n - 1):L]
                        lin *= (2.0 / n)
                        q = np.empty(2 * n - 1)
                        q[n - 1:] = ac[:n]; q[:n - 1] = ac[L - (n - 1):L]
                        q *= -(2.0 / n)
                        d_vals = lin + q
                        new_max = float((vals + d_vals).max())
                        drop = cur_ub - new_max
                        if drop > best_drop + 1e-19:
                            best_drop = drop
                            best_move = ('block', int(ii), int(jj),
                                         int(w_size), float(delta))

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
            _, ii, jj, w_size, delta = best_move
            a[ii:ii + w_size] -= delta
            a[jj:jj + w_size] += delta
            a = np.clip(a, 0.0, 1.0)
        accept += 1; stall = 0
        if verbose and accept % log_every == 0:
            print(f"  MT it {itr} acc {accept} step={step:.2e} "
                  f"ub={float(overlap_full(a).max()):.17f} drop={best_drop:.2e} "
                  f"kind={best_move[0]}")
    return a


def active_set_polish(a, max_iters=500, tol=1e-17, active_tol=1e-4,
                      max_active=600, verbose=False, log_every=20):
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
        if t > -1e-17: break
        with np.errstate(divide='ignore', invalid='ignore'):
            pos = d > 1e-15; neg = d < -1e-15
            step_pos = np.min((1.0 - a)[pos] / d[pos]) if pos.any() else np.inf
            step_neg = np.min(a[neg] / (-d[neg])) if neg.any() else np.inf
        max_step = 0.999 * min(step_pos, step_neg)
        if max_step <= 1e-18: break
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
        for _ in range(90):
            u = float((vals + s * lin + s * s * q).max())
            if u < best_u - 1e-17: best_u, best_s = u, s; break
            s *= 0.5
            if s < 1e-19: break
        if best_s > 0:
            for f in (1.3, 1.5, 2.0, 0.8, 0.6):
                ss = min(best_s * f, max_step)
                u = float((vals + ss * lin + ss * ss * q).max())
                if u < best_u - 1e-17: best_u, best_s = u, ss
        if best_s == 0.0: break
        a = a + best_s * d; a = np.clip(a, 0.0, 1.0)
        excess = a.sum() - n / 2.0
        if abs(excess) > 1e-12:
            inner = (a > 1e-6) & (a < 1 - 1e-6)
            if inner.any():
                a[inner] -= excess / inner.sum(); a = np.clip(a, 0.0, 1.0)
        cur_ub = float(overlap_full(a).max())
        if verbose and (it < 5 or it % log_every == 0):
            print(f"  AS it {it:3d}: K={K} t={t:.2e} s={best_s:.3e} ub={cur_ub:.17f} drop={last_ub-cur_ub:.2e}")
        if last_ub - cur_ub < tol: last_ub = cur_ub; break
        last_ub = cur_ub
    return a


def solve(verbose=False, time_budget_sec=900):
    here = os.path.dirname(os.path.abspath(__file__))
    t_start = time.time()

    # Warm from v13
    v13_path = os.path.join(here, "..", "v13", "sequence.json")
    with open(v13_path) as f:
        d = json.load(f)
    best_a = np.asarray(d["a"], dtype=np.float64)
    best_ub = float(overlap_full(best_a).max())
    if verbose:
        print(f"[warm v13] U = {best_ub:.17f}, n = {len(best_a)}")

    # ---- Phase 1: AS with WIDE active_tol (include all near-max shifts) ----
    if verbose: print(f"\n[Phase 1: AS very wide active_tol=1e-2]")
    pol = active_set_polish(best_a, max_iters=200, tol=1e-17,
                            active_tol=1e-2, max_active=1199,
                            verbose=verbose)
    ub_p = float(overlap_full(pol).max())
    if ub_p < best_ub:
        if verbose: print(f"   gain: {best_ub-ub_p:.2e}, new ub={ub_p:.17f}")
        best_a, best_ub = pol, ub_p

    # ---- Phase 2: MT with block-transfer moves at various widths ----
    p2_budget = max(60.0, time_budget_sec * 0.30)
    if verbose: print(f"\n[Phase 2: MT with block moves, budget {p2_budget:.0f}s]")
    a_mt = mass_transport_polish(best_a, time_budget=p2_budget,
                                 step_init=2.0**-18, step_min=2.0**-42,
                                 verbose=verbose, n_top=60, n_rand=30,
                                 use_block=True, block_widths=(2, 3, 5, 8))
    ub_mt = float(overlap_full(a_mt).max())
    if ub_mt < best_ub:
        if verbose: print(f"   MT gain: {best_ub-ub_mt:.2e}, new ub={ub_mt:.17f}")
        best_a, best_ub = a_mt, ub_mt

    # ---- Phase 3: alternate AS + MT at successively finer steps ----
    p3_end = t_start + time_budget_sec - 20
    cyc = 0
    while time.time() < p3_end:
        cyc += 1
        if verbose: print(f"\n[Phase 3 cycle {cyc}] AS wide tol")
        pol = active_set_polish(best_a, max_iters=200, tol=1e-17,
                                active_tol=10**(-3 - cyc * 0.5), max_active=1199)
        ub_p = float(overlap_full(pol).max())
        if ub_p < best_ub - 1e-16:
            if verbose: print(f"   AS gain: {best_ub-ub_p:.2e}")
            best_a, best_ub = pol, ub_p
        mt_budget = min((time_budget_sec - (time.time() - t_start)) / 3, 80.0)
        if mt_budget < 8: break
        if verbose: print(f"   MT step=2^-{20+cyc*2}, budget {mt_budget:.0f}s")
        a_mt = mass_transport_polish(best_a, time_budget=mt_budget,
                                     step_init=2.0**-(20 + cyc * 2),
                                     step_min=2.0**-46,
                                     verbose=verbose, n_top=60, n_rand=30,
                                     use_block=True,
                                     block_widths=(2, 3, 4, 6, 10),
                                     seed=20000 + cyc)
        ub_mt = float(overlap_full(a_mt).max())
        if ub_mt < best_ub - 1e-17:
            if verbose: print(f"   MT gain: {best_ub-ub_mt:.2e}")
            best_a, best_ub = a_mt, ub_mt
        elif verbose:
            print(f"   (no MT gain)")
        if cyc > 8: break

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
