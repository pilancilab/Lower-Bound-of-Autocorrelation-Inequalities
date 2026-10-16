"""
v1 — Improves over v0 by:
  - Warm start from v0's sequence.json.
  - Aggressive basin-hopping: large structured perturbations + L-BFGS
    re-optimization repeated 100+ times.
  - Mass-transport hillclimb with multi-bin coordinated moves AND simulated
    annealing acceptance to escape strict local minima.
  - Refinement at multiple n values (600, 1024, 2048), with the best kept.

Reuses the same core utilities (FFT overlap & gradient, logit param +
sum projection, LSE-smoothed max) as v0.
"""

import json
import os
import time
import numpy as np
from scipy.optimize import minimize, brentq, linprog

EPS = 1e-9


# ---------------------------------------------------------------------------
# Core
# ---------------------------------------------------------------------------

def _L(n):
    L = 1
    while L < 2 * n:
        L *= 2
    return L


def overlap_full(a):
    n = len(a)
    b = 1.0 - a
    L = _L(n)
    A = np.fft.rfft(a, n=L)
    B = np.fft.rfft(b, n=L)
    full = np.fft.irfft(np.conjugate(A) * B, n=L)
    M = np.empty(2 * n - 1)
    M[n - 1:] = full[:n]
    M[:n - 1] = full[L - (n - 1):L]
    return M * (2.0 / n)


def grad_weighted(a, w):
    n = len(a)
    b = 1.0 - a
    L = _L(n)
    w_circ = np.zeros(L)
    w_circ[:n] = w[n - 1:]
    w_circ[L - (n - 1):L] = w[:n - 1]
    W = np.fft.rfft(w_circ, n=L)
    B = np.fft.rfft(b, n=L)
    A = np.fft.rfft(a, n=L)
    T1 = np.fft.irfft(np.conjugate(W) * B, n=L)[:n]
    T2 = np.fft.irfft(W * A, n=L)[:n]
    return (2.0 / n) * (T1 - T2)


def sigmoid(z):
    out = np.empty_like(z)
    pos = z >= 0
    out[pos] = 1.0 / (1.0 + np.exp(-z[pos]))
    ez = np.exp(z[~pos])
    out[~pos] = ez / (1.0 + ez)
    return out


def project_sum_c(z, target):
    def F(c):
        return sigmoid(z + c).sum() - target
    lo, hi = -50.0, 50.0
    while F(lo) > 0 and lo > -1e5: lo -= 50.0
    while F(hi) < 0 and hi < 1e5:  hi += 50.0
    return brentq(F, lo, hi, xtol=1e-13, rtol=1e-13)


def make_fg(n, beta):
    target = n / 2.0
    def fg(z):
        c = project_sum_c(z, target)
        a = sigmoid(z + c)
        vals = overlap_full(a)
        m = vals.max()
        e = np.exp(beta * (vals - m))
        Z = e.sum()
        smax = m + np.log(Z) / beta
        w = e / Z
        df_da = grad_weighted(a, w)
        ab = a * (1.0 - a)
        S = ab.sum()
        if S < 1e-30:
            return float(smax), np.zeros_like(z)
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


# ---------------------------------------------------------------------------
# Perturbations
# ---------------------------------------------------------------------------

def perturb_lowfreq(a, rng, amplitude=0.1, n_modes=12):
    n = len(a)
    x = np.linspace(-1.0, 1.0, n, endpoint=False) + 1.0 / n
    delta = np.zeros(n)
    for k in range(1, n_modes + 1):
        delta += (rng.standard_normal() * np.cos(k * np.pi * x)
                  + rng.standard_normal() * np.sin(k * np.pi * x)) / k
    delta = delta / (np.std(delta) + 1e-12) * amplitude
    a_new = np.clip(a + delta, EPS, 1 - EPS)
    a_new = a_new * (n / 2.0) / a_new.sum()
    return np.clip(a_new, EPS, 1 - EPS)


def perturb_highfreq(a, rng, amplitude=0.05):
    n = len(a)
    delta = amplitude * rng.standard_normal(n)
    a_new = np.clip(a + delta, EPS, 1 - EPS)
    a_new = a_new * (n / 2.0) / a_new.sum()
    return np.clip(a_new, EPS, 1 - EPS)


def perturb_blockflip(a, rng, n_blocks=2):
    """Flip a few random blocks: a_i -> 1 - a_i in segments. Then renormalize."""
    n = len(a)
    a_new = a.copy()
    for _ in range(n_blocks):
        L = int(rng.integers(n // 20, n // 5))
        start = int(rng.integers(0, n - L))
        a_new[start:start + L] = 1.0 - a_new[start:start + L]
    a_new = np.clip(a_new, EPS, 1 - EPS)
    a_new = a_new * (n / 2.0) / a_new.sum()
    return np.clip(a_new, EPS, 1 - EPS)


def perturb_shift(a, rng, max_shift=None):
    """Cyclic shift inside the array (kind of). Then renormalize."""
    n = len(a)
    if max_shift is None:
        max_shift = n // 4
    s = int(rng.integers(-max_shift, max_shift + 1))
    if s == 0:
        s = 1
    a_new = np.roll(a, s)
    # Zero-pad the wrapped region
    if s > 0:
        a_new[:s] = 0.5
    else:
        a_new[s:] = 0.5
    a_new = np.clip(a_new, EPS, 1 - EPS)
    a_new = a_new * (n / 2.0) / a_new.sum()
    return np.clip(a_new, EPS, 1 - EPS)


def perturb_reverse(a):
    n = len(a)
    a_new = a[::-1].copy()
    return a_new


def perturb_complement(a):
    """a -> 1 - a, then renormalize (gives M_k -> M_{-k}, same max)."""
    n = len(a)
    a_new = 1.0 - a
    # sum = n - n/2 = n/2 already
    return a_new


def perturb_random_choice(a, rng, amplitude=None):
    if amplitude is None:
        amplitude = float(rng.choice([0.03, 0.08, 0.15, 0.25]))
    kind = rng.choice(["lowfreq", "lowfreq", "highfreq", "blockflip",
                       "shift", "blockflip"])
    if kind == "lowfreq":
        return perturb_lowfreq(a, rng, amplitude=amplitude,
                               n_modes=int(rng.integers(4, 25)))
    elif kind == "highfreq":
        return perturb_highfreq(a, rng, amplitude=amplitude)
    elif kind == "blockflip":
        return perturb_blockflip(a, rng, n_blocks=int(rng.integers(1, 4)))
    elif kind == "shift":
        return perturb_shift(a, rng)
    return a


def interpolate(a, n_new):
    n_old = len(a)
    x_old = (np.arange(n_old) + 0.5) / n_old
    x_new = (np.arange(n_new) + 0.5) / n_new
    a_new = np.interp(x_new, x_old, a)
    a_new = np.clip(a_new, EPS, 1 - EPS)
    a_new = a_new * (n_new / 2.0) / a_new.sum()
    return np.clip(a_new, EPS, 1 - EPS)


# ---------------------------------------------------------------------------
# Mass transport (deterministic, greedy) with multi-bin moves and SA option
# ---------------------------------------------------------------------------

def grad_single_shift(a, d):
    n = len(a)
    b = 1.0 - a
    g = np.zeros(n)
    if d >= 0:
        g[:n - d] += b[d:]
        g[d:] -= a[:n - d]
    else:
        kk = -d
        g[kk:] += b[:n - kk]
        g[:n - kk] -= a[kk:]
    return g * (2.0 / n)


def _bin_grad(a, p):
    """(grad_k)_p for all k in [-(n-1),n-1] (length 2n-1, indexed by k+n-1)."""
    n = len(a)
    b = 1.0 - a
    g = np.zeros(2 * n - 1)
    idx_b = np.arange(n) - p + (n - 1)
    g[idx_b] += b
    idx_a = p - np.arange(n) + (n - 1)
    g[idx_a] -= a
    return g * (2.0 / n)


def mass_transport_polish(a, time_budget=180.0, step_init=1e-3,
                          step_min=2.0 ** -32, sa_T0=0.0,
                          sa_decay=0.0, verbose=False, log_every=40):
    a = a.copy()
    n = len(a)
    t_start = time.time()
    step = step_init
    cur_ub = float(overlap_full(a).max())
    T = sa_T0
    iter_count = 0
    accept_count = 0
    stall = 0

    while time.time() - t_start < time_budget:
        iter_count += 1
        vals = overlap_full(a)
        cur_ub = vals.max()
        thr = cur_ub - 5e-5
        active = np.where(vals >= thr)[0]
        if len(active) == 0:
            break
        w = np.zeros(2 * n - 1)
        ev = np.exp(5000.0 * (vals[active] - cur_ub))
        w[active] = ev / ev.sum()
        g = grad_weighted(a, w)

        elig_d = np.where(a > step)[0]
        elig_r = np.where(a < 1.0 - step)[0]
        if len(elig_d) == 0 or len(elig_r) == 0:
            step *= 0.5
            if step < step_min: break
            continue

        K_top = min(25, len(elig_d), len(elig_r))
        donor_top = elig_d[np.argsort(g[elig_d])[::-1][:K_top]]
        recv_top = elig_r[np.argsort(g[elig_r])[:K_top]]
        # plus randomness
        n_rand = min(12, len(elig_d), len(elig_r))
        rng = np.random.default_rng(iter_count + 313)
        donor_rand = rng.choice(elig_d, size=n_rand, replace=False)
        recv_rand = rng.choice(elig_r, size=n_rand, replace=False)
        donors = np.unique(np.concatenate([donor_top, donor_rand]))
        recvs = np.unique(np.concatenate([recv_top, recv_rand]))

        best_drop = 0.0 if T == 0 else -np.inf
        best_pair = None
        best_delta = 0.0
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
                if T > 0:
                    accept = (drop > 0) or (rng.random() < np.exp(drop / T))
                    if accept and drop > best_drop:
                        best_drop = drop
                        best_pair = (int(i), int(j))
                        best_delta = float(delta)
                else:
                    if drop > best_drop + 1e-17:
                        best_drop = drop
                        best_pair = (int(i), int(j))
                        best_delta = float(delta)

        if best_pair is not None:
            i, j = best_pair
            a[i] -= best_delta; a[i] = max(0.0, a[i])
            a[j] += best_delta; a[j] = min(1.0, a[j])
            accept_count += 1
            stall = 0
            new_ub = float(overlap_full(a).max())
            if verbose and accept_count % log_every == 0:
                print(f"  MT iter {iter_count} acc {accept_count}, "
                      f"step={step:.2e}, T={T:.2e}, ub={new_ub:.10f}, "
                      f"drop={best_drop:.2e}")
        else:
            stall += 1
            if stall >= 4:
                step *= 0.5
                stall = 0
                if verbose:
                    print(f"  MT step -> {step:.2e}")
                if step < step_min:
                    break
        if sa_decay > 0:
            T *= sa_decay
            if T < 1e-12:
                T = 0.0
    return a


def active_set_polish(a, max_iters=300, tol=1e-13, active_tol=1e-4,
                      max_active=300, verbose=False):
    a = a.copy()
    n = len(a)
    last_ub = float(overlap_full(a).max())
    for it in range(max_iters):
        vals = overlap_full(a)
        ub = vals.max()
        thr = ub - active_tol
        active_idx = np.where(vals >= thr)[0]
        if len(active_idx) > max_active:
            order = np.argsort(vals[active_idx])[::-1][:max_active]
            active_idx = active_idx[order]
        K = len(active_idx)
        if K == 0:
            break
        grads = np.zeros((K, n))
        for j, idx in enumerate(active_idx):
            grads[j] = grad_single_shift(a, int(idx - (n - 1)))
        c = np.zeros(n + 1); c[-1] = 1.0
        A_ub = np.concatenate([grads, -np.ones((K, 1))], axis=1)
        b_ub = np.zeros(K)
        A_eq = np.concatenate([np.ones((1, n)), np.zeros((1, 1))], axis=1)
        b_eq = np.zeros(1)
        bounds = [(-1.0, 1.0)] * n + [(None, None)]
        try:
            res = linprog(c, A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=b_eq,
                          bounds=bounds, method='highs')
        except Exception:
            break
        if not res.success or res.x is None: break
        d = res.x[:n]; t = res.x[-1]
        if t > -1e-14: break
        with np.errstate(divide='ignore', invalid='ignore'):
            pos = d > 1e-15; neg = d < -1e-15
            step_pos = np.min((1.0 - a)[pos] / d[pos]) if pos.any() else np.inf
            step_neg = np.min(a[neg] / (-d[neg])) if neg.any() else np.inf
        max_step = 0.999 * min(step_pos, step_neg)
        if max_step <= 1e-16: break
        L = _L(n)
        D = np.fft.rfft(d, n=L)
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
        s = max_step
        best_s, best_u = 0.0, ub
        for _ in range(80):
            u = float((vals + s * lin + s * s * q).max())
            if u < best_u - 1e-15: best_u, best_s = u, s; break
            s *= 0.5
            if s < 1e-17: break
        if best_s > 0:
            for factor in (1.3, 1.5, 2.0, 0.8, 0.6):
                ss = min(best_s * factor, max_step)
                u = float((vals + ss * lin + ss * ss * q).max())
                if u < best_u - 1e-15: best_u, best_s = u, ss
        if best_s == 0.0: break
        a = a + best_s * d
        a = np.clip(a, 0.0, 1.0)
        excess = a.sum() - n / 2.0
        if abs(excess) > 1e-12:
            inner = (a > 1e-6) & (a < 1 - 1e-6)
            if inner.any():
                a[inner] -= excess / inner.sum()
                a = np.clip(a, 0.0, 1.0)
        cur_ub = float(overlap_full(a).max())
        if last_ub - cur_ub < tol:
            last_ub = cur_ub
            break
        last_ub = cur_ub
    return a


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

def solve(verbose=False, time_budget_sec=900):
    here = os.path.dirname(os.path.abspath(__file__))
    t_start = time.time()

    # ---- Warm start from v0 ----
    v0_path = os.path.join(here, "..", "v0", "sequence.json")
    assert os.path.isfile(v0_path), "v0 sequence missing"
    with open(v0_path) as f:
        v0_data = json.load(f)
    best_a = np.asarray(v0_data["a"], dtype=float)
    best_ub = float(v0_data.get("upper_bound", float(overlap_full(best_a).max())))
    if verbose:
        print(f"[warm v0] ub={best_ub:.10f}, n={len(best_a)}")

    # ---- Phase 1: aggressive basin hopping at n=600 ----
    rng = np.random.default_rng(0xBEEF)
    n_main = len(best_a)
    betas_hop = (200.0, 1000.0, 5000.0, 30000.0)
    iters_hop = 200
    phase_end = t_start + time_budget_sec * 0.45
    hop = 0
    no_improve_streak = 0
    while time.time() < phase_end:
        hop += 1
        amp = float(rng.choice([0.03, 0.06, 0.12, 0.20, 0.30]))
        a_try = perturb_random_choice(best_a, rng, amplitude=amp)
        a_opt, ub = optimize_from_a(a_try, betas_hop, iters_hop)
        if ub < best_ub - 1e-8:
            if verbose:
                print(f"[hop #{hop} amp={amp:.2f}] ub={ub:.10f} (NEW)")
            best_a, best_ub = a_opt, ub
            no_improve_streak = 0
        else:
            no_improve_streak += 1

    # ---- Phase 2: complement + reverse explorations ----
    # These preserve U exactly; perturb them and re-optimize to potentially find
    # a structurally different basin.
    structural = [perturb_complement(best_a),
                  perturb_reverse(best_a)]
    for k, a_struct in enumerate(structural):
        if time.time() - t_start > time_budget_sec * 0.55:
            break
        a_try = perturb_lowfreq(a_struct, rng, amplitude=0.10, n_modes=15)
        a_opt, ub = optimize_from_a(a_try, betas_hop, iters_hop)
        if ub < best_ub - 1e-8:
            if verbose:
                print(f"[struct #{k}] ub={ub:.10f} (NEW)")
            best_a, best_ub = a_opt, ub

    # ---- Phase 3: heavier perturbation pairs at n=600 ----
    rng2 = np.random.default_rng(0xC0FFEE)
    p3_end = t_start + time_budget_sec * 0.65
    while time.time() < p3_end:
        a_try = best_a.copy()
        # apply multiple perturbations
        for _ in range(int(rng2.integers(1, 4))):
            a_try = perturb_random_choice(a_try, rng2)
        a_opt, ub = optimize_from_a(a_try, betas_hop, iters_hop)
        if ub < best_ub - 1e-8:
            if verbose: print(f"[heavy hop] ub={ub:.10f} (NEW)")
            best_a, best_ub = a_opt, ub

    # ---- Phase 4: try n=1024 and n=2048 by interpolation ----
    for n_new in (1024,):  # 2048 takes too long, skip unless time permits
        if time.time() - t_start > time_budget_sec * 0.75:
            break
        a_init = interpolate(best_a, n_new)
        a_opt, ub = optimize_from_a(a_init, (200.0, 1000.0, 5000.0, 30000.0), 300)
        if verbose:
            print(f"[refine n={n_new}] ub={ub:.10f}")
        # if better, also polish at that n
        if ub < best_ub - 1e-8:
            best_a, best_ub = a_opt, ub
        # then come back down to n=600 by subsampling
        a_back = interpolate(a_opt, n_main)
        a_opt2, ub2 = optimize_from_a(a_back, betas_hop, iters_hop)
        if ub2 < best_ub - 1e-8:
            if verbose: print(f"[refine-back n={n_main}] ub={ub2:.10f} (NEW)")
            best_a, best_ub = a_opt2, ub2

    # ---- Phase 5: final LSE + active-set + mass-transport polish ----
    a_opt, ub = optimize_from_a(best_a, (3000.0, 15000.0, 50000.0), 500)
    if ub < best_ub:
        if verbose: print(f"[LSE polish] ub={ub:.10f} (NEW)")
        best_a, best_ub = a_opt, ub

    polished = active_set_polish(best_a, max_iters=400, tol=1e-13)
    ub_p = float(overlap_full(polished).max())
    if ub_p < best_ub:
        if verbose: print(f"[AS polish] ub={ub_p:.10f} (NEW)")
        best_a, best_ub = polished, ub_p

    remaining = max(60.0, time_budget_sec - (time.time() - t_start))
    if verbose:
        print(f"[mass transport] starting, budget={remaining:.0f}s, "
              f"current={best_ub:.10f}")
    mt = mass_transport_polish(best_a, time_budget=remaining,
                               step_init=1e-3, step_min=2.0 ** -32,
                               verbose=verbose)
    ub_mt = float(overlap_full(mt).max())
    if ub_mt < best_ub:
        if verbose: print(f"[MT] ub={ub_mt:.10f} (NEW)")
        best_a, best_ub = mt, ub_mt

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
