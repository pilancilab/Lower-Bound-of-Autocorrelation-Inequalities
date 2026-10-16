"""
v6 — Explores higher-n refinement of v5 (warm start), and asymmetric inits.

We interpolate v5's n=600 solution to n=1200 (exact 2x repeat) and to n=900
(linear interp) to see if more degrees of freedom help the mass-transport
polish escape its current local minimum. Then asymmetric inits via random
shifts of the bump pattern.
"""

import json, os, time
import numpy as np
from scipy.optimize import minimize, brentq, linprog

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


def normalize(a):
    n = len(a)
    a = np.clip(a, EPS, 1 - EPS)
    a = a * (n / 2.0) / a.sum()
    return np.clip(a, EPS, 1 - EPS)


def upsample_repeat(a, factor):
    """Repeat each bin 'factor' times, preserving the step function shape exactly."""
    return np.repeat(a, factor)


def downsample_mean(a, factor):
    n = len(a)
    assert n % factor == 0
    return a.reshape(n // factor, factor).mean(axis=1)


def perturb_lowfreq(a, rng, amplitude=0.1, n_modes=12):
    n = len(a)
    x = np.linspace(-1.0, 1.0, n, endpoint=False) + 1.0 / n
    delta = np.zeros(n)
    for k in range(1, n_modes + 1):
        delta += (rng.standard_normal() * np.cos(k * np.pi * x)
                  + rng.standard_normal() * np.sin(k * np.pi * x)) / k
    delta = delta / (np.std(delta) + 1e-12) * amplitude
    return normalize(a + delta)


def active_set_polish(a, max_iters=400, tol=1e-14, active_tol=1e-4,
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


def mass_transport_polish(a, time_budget=180.0, step_init=1e-3,
                          step_min=2.0 ** -36, verbose=False, log_every=120):
    a = a.copy(); n = len(a); t_start = time.time()
    step = step_init; accept = 0; stall = 0; itr = 0
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
        K_top = min(40, len(elig_d), len(elig_r))
        donor_top = elig_d[np.argsort(g[elig_d])[::-1][:K_top]]
        recv_top = elig_r[np.argsort(g[elig_r])[:K_top]]
        rng = np.random.default_rng(itr + 313)
        donor_rand = rng.choice(elig_d, size=min(20, len(elig_d)), replace=False)
        recv_rand = rng.choice(elig_r, size=min(20, len(elig_r)), replace=False)
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
                if drop > best_drop + 1e-17:
                    best_drop = drop
                    best_pair = (int(i), int(j))
                    best_delta = float(delta)
        if best_pair is not None:
            i, j = best_pair
            a[i] -= best_delta; a[i] = max(0.0, a[i])
            a[j] += best_delta; a[j] = min(1.0, a[j])
            accept += 1; stall = 0
            if verbose and accept % log_every == 0:
                print(f"  MT it {itr} acc {accept} step={step:.2e} "
                      f"ub={float(overlap_full(a).max()):.17f} drop={best_drop:.2e}")
        else:
            stall += 1
            if stall >= 4:
                step *= 0.5; stall = 0
                if step < step_min: break
    return a


# ---- driver ----------------------------------------------------------------

def solve(verbose=False, time_budget_sec=900):
    here = os.path.dirname(os.path.abspath(__file__))
    t_start = time.time()
    betas_short = (300.0, 1500.0, 8000.0, 40000.0)
    iters_short = 200

    v5_path = os.path.join(here, "..", "v5", "sequence.json")
    assert os.path.isfile(v5_path)
    with open(v5_path) as f:
        v5_data = json.load(f)
    a_v5 = np.asarray(v5_data["a"], dtype=float)
    ub_v5 = float(v5_data.get("upper_bound", float(overlap_full(a_v5).max())))
    if verbose:
        print(f"[warm v5] ub={ub_v5:.17f}, n={len(a_v5)}")

    best_a, best_ub = a_v5, ub_v5

    # ---- Phase 1: try n=1200 (exact 2x upsample, same step function) ----
    if time.time() - t_start < time_budget_sec * 0.10:
        a_up = upsample_repeat(a_v5, 2)
        ub_up = float(overlap_full(a_up).max())
        if verbose:
            print(f"[upsample n=1200] ub={ub_up:.17f} (should equal v5)")
        # The discrete max at n=1200 may differ slightly from n=600 because
        # the additional shifts at n=1200 may probe new x values.
        # Polish at n=1200
        a_p, ub_p = optimize_from_a(a_up, betas_short, iters_short)
        if verbose: print(f"[n=1200 LSE] {ub_p:.17f}")
        if ub_p < best_ub - 1e-9:
            best_a, best_ub = a_p, ub_p
            if verbose: print(f"   (NEW BEST at n=1200)")
        # Polish with MT
        a_p = active_set_polish(a_p, max_iters=300, tol=1e-14)
        ub_p2 = float(overlap_full(a_p).max())
        if ub_p2 < best_ub: best_a, best_ub = a_p, ub_p2
        # MT at n=1200
        mt_budget = time_budget_sec * 0.25
        if verbose: print(f"[n=1200 MT] starting, budget {mt_budget:.0f}s")
        a_mt = mass_transport_polish(a_p, time_budget=mt_budget,
                                     step_init=1e-4, step_min=2.0 ** -36,
                                     verbose=verbose)
        ub_mt = float(overlap_full(a_mt).max())
        if verbose: print(f"[n=1200 MT done] {ub_mt:.17f}")
        if ub_mt < best_ub:
            best_a, best_ub = a_mt, ub_mt
            if verbose: print(f"   (NEW BEST at n=1200)")

    # ---- Phase 2: try n=900 ----
    if time.time() - t_start < time_budget_sec * 0.50:
        x_old = (np.arange(len(a_v5)) + 0.5) / len(a_v5)
        x_new = (np.arange(900) + 0.5) / 900
        a_900 = np.interp(x_new, x_old, a_v5)
        a_900 = normalize(a_900)
        a_p, ub_p = optimize_from_a(a_900, betas_short, iters_short)
        if verbose: print(f"[n=900 LSE] {ub_p:.17f}")
        if ub_p < best_ub - 1e-9:
            best_a, best_ub = a_p, ub_p
        a_p = active_set_polish(a_p, max_iters=200, tol=1e-14)
        ub_p2 = float(overlap_full(a_p).max())
        if ub_p2 < best_ub: best_a, best_ub = a_p, ub_p2

    # ---- Phase 3: heavy basin hopping at current best n ----
    p3_end = t_start + time_budget_sec * 0.75
    rng = np.random.default_rng(0xCAFE)
    hop = 0
    n_cur = len(best_a)
    while time.time() < p3_end:
        hop += 1
        amp = float(rng.choice([0.05, 0.15, 0.30, 0.45]))
        modes = int(rng.integers(3, 25))
        a_try = perturb_lowfreq(best_a, rng, amplitude=amp, n_modes=modes)
        a_opt, ub = optimize_from_a(a_try, betas_short, iters_short)
        if ub < best_ub - 1e-9:
            if verbose: print(f"[hop #{hop} amp={amp:.2f}] {ub:.17f} (NEW)")
            best_a, best_ub = a_opt, ub

    # ---- Phase 4: final polish ----
    a_p, ub_p = optimize_from_a(best_a, (5000.0, 25000.0, 100000.0), 600)
    if ub_p < best_ub:
        if verbose: print(f"[final LSE] {ub_p:.17f} (NEW)")
        best_a, best_ub = a_p, ub_p
    p2 = active_set_polish(best_a, max_iters=400, tol=1e-14)
    ub_p2 = float(overlap_full(p2).max())
    if ub_p2 < best_ub:
        if verbose: print(f"[final AS] {ub_p2:.17f} (NEW)")
        best_a, best_ub = p2, ub_p2
    remaining = max(40.0, time_budget_sec - (time.time() - t_start))
    if verbose: print(f"[final MT] budget {remaining:.0f}s")
    a_mt = mass_transport_polish(best_a, time_budget=remaining,
                                 step_init=1e-5, step_min=2.0 ** -36,
                                 verbose=verbose)
    ub_mt = float(overlap_full(a_mt).max())
    if ub_mt < best_ub:
        if verbose: print(f"[final MT] {ub_mt:.17f} (NEW)")
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
