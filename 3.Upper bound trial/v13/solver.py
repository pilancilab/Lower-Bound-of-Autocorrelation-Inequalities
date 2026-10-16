"""
v13 — Warm start from the actual online EinsteinArena SOTA:
DarwinAgent8427's sequence (thread #86) at U = 0.38087031046844222.

This is strictly better than:
  - Together AI 2026 github file (U = 0.38087031058621990) by 1.18e-10
  - My v12 result            (U = 0.38087031050250464) by 3.40e-11

Pipeline: AS LP polish + long dyadic MT polish at fine step. Multiple
cycles to alternately re-equalize and squeeze.
"""

import json, os, time, sys
import numpy as np
from scipy.optimize import minimize, brentq, linprog

# Make sota_ref importable
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))

EPS = 1e-9


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


def mass_transport_polish(a, time_budget=180.0, step_init=2.0**-20,
                          step_min=2.0**-44, verbose=False, log_every=200,
                          n_top=50, n_rand=25, seed=777):
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
        best_drop = 0.0; best_pair = None; best_delta = 0.0
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
                    best_pair = (int(i), int(j)); best_delta = float(delta)
        if best_pair is not None:
            i, j = best_pair
            a[i] -= best_delta; a[i] = max(0.0, a[i])
            a[j] += best_delta; a[j] = min(1.0, a[j])
            accept += 1; stall = 0
            if verbose and accept % log_every == 0:
                print(f"  MT it {itr} acc {accept} step={step:.2e} "
                      f"ub={float(overlap_full(a).max()):.17f}")
        else:
            stall += 1
            if stall >= 4:
                step *= 0.5; stall = 0
                if step < step_min: break
    return a


def solve(verbose=False, time_budget_sec=900):
    here = os.path.dirname(os.path.abspath(__file__))
    t_start = time.time()

    # Load DarwinAgent8427 sequence
    src = os.path.join(here, "..", "sota_ref", "darwin_agent_8427.json")
    with open(src) as f:
        d = json.load(f)
    best_a = np.asarray(d["a"], dtype=np.float64)
    best_ub = float(overlap_full(best_a).max())
    if verbose:
        print(f"[warm DarwinAgent8427] U = {best_ub:.17f}, n = {len(best_a)}")
        print(f"   (claimed by author: {d['claimed_score']})")

    # Repeated AS+MT cycles, tightening tol and step each cycle
    cycles = [
        # (active_tol, mt_step_init, mt_budget_frac)
        (1e-4, 2.0**-18, 0.20),
        (1e-5, 2.0**-22, 0.15),
        (1e-6, 2.0**-26, 0.15),
        (1e-7, 2.0**-30, 0.20),
        (1e-8, 2.0**-34, 0.30),
    ]
    for cyc, (atol, mt_step, mt_frac) in enumerate(cycles):
        if time.time() - t_start > time_budget_sec - 20: break
        if verbose: print(f"\n[Cycle {cyc}: AS atol={atol:.0e}] cur {best_ub:.17f}")
        pol = active_set_polish(best_a, max_iters=300, tol=1e-17,
                                active_tol=atol, max_active=600,
                                verbose=verbose)
        ub_p = float(overlap_full(pol).max())
        if ub_p < best_ub:
            if verbose: print(f"   AS gain: {best_ub - ub_p:.2e}, ub={ub_p:.17f}")
            best_a, best_ub = pol, ub_p
        mt_budget = min(time_budget_sec * mt_frac,
                        time_budget_sec - (time.time() - t_start) - 5)
        if mt_budget < 10: continue
        if verbose: print(f"[Cycle {cyc}: MT step={mt_step:.2e}, budget {mt_budget:.0f}s]")
        a_mt = mass_transport_polish(best_a, time_budget=mt_budget,
                                     step_init=mt_step, step_min=2.0**-44,
                                     verbose=verbose, n_top=60, n_rand=30,
                                     seed=10001 + cyc)
        ub_mt = float(overlap_full(a_mt).max())
        if ub_mt < best_ub:
            if verbose: print(f"   MT gain: {best_ub - ub_mt:.2e}, ub={ub_mt:.17f}")
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
