"""
v0 — constructive upper bound for the Erdős minimum overlap problem.

Pipeline (per proposal.md):
  1. Multi-init L-BFGS on LSE-smoothed max with logit + sum-projection.
  2. Multi-resolution refinement (lower n -> higher n by interpolation).
  3. Active-set sequential LP for descent direction on tied shifts +
     exact quadratic line search.
  4. Final mass-transport hillclimb on a dyadic grid (sum- and box-preserving
     by construction), inspired by the current SOTA approach
     (TuringAgent9811, 0.38087031 on n=600).

The script defines solve() returning a dict per spec.
"""

import json
import os
import time
import numpy as np
from scipy.optimize import minimize, brentq, linprog

EPS = 1e-9


# ---------------------------------------------------------------------------
# Core FFT-based overlap and gradient
# ---------------------------------------------------------------------------

def _L(n):
    L = 1
    while L < 2 * n:
        L *= 2
    return L


def overlap_full(a):
    """vals[d + n - 1] = M_d = (2/n) sum_i a_i (1 - a_{i+d}),  d in [-(n-1),n-1]."""
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
    """Gradient of sum_d w[d+n-1] * M_d w.r.t. a (FFT, O(n log n))."""
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


def grad_single_shift(a, d):
    """∂M_d/∂a."""
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


# ---------------------------------------------------------------------------
# Logit parameterization with sum-projection
# ---------------------------------------------------------------------------

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
    z = np.log(a) - np.log(1 - a)
    z = z + project_sum_c(z, n / 2.0)
    z = run_lbfgs(z, n, betas, iters)
    a_out = sigmoid(z + project_sum_c(z, n / 2.0))
    return a_out, float(overlap_full(a_out).max())


# ---------------------------------------------------------------------------
# Initializations and perturbations
# ---------------------------------------------------------------------------

def make_init(n, kind, rng):
    x = np.linspace(-1.0, 1.0, n, endpoint=False) + 1.0 / n
    if kind == "cosine":
        a = 0.5 + 0.4 * np.cos(np.pi * x)
    elif kind == "double_cos":
        a = 0.5 + 0.3 * np.cos(2 * np.pi * x)
    elif kind == "triple_cos":
        a = 0.5 + 0.2 * np.cos(np.pi * x) + 0.2 * np.cos(3 * np.pi * x)
    elif kind == "tent":
        a = 1.0 - np.abs(x)
    elif kind == "bell":
        a = np.exp(-2.0 * x * x)
    elif kind == "bell_narrow":
        a = np.exp(-5.0 * x * x)
    elif kind == "bell_wide":
        a = np.exp(-0.5 * x * x)
    elif kind == "raised_cos":
        a = 0.6 - 0.5 * np.cos(np.pi * x) ** 2
    elif kind == "two_block":
        a = np.zeros(n)
        a[n // 10: 4 * n // 10] = 1.0
        a[6 * n // 10: 9 * n // 10] = 1.0
        a += 0.03 * rng.standard_normal(n)
    elif kind == "rand_smooth":
        z = rng.standard_normal(n)
        k = np.fft.fftfreq(n) * n
        mask = (np.abs(k) <= max(3, n // 25)).astype(float)
        z = np.real(np.fft.ifft(np.fft.fft(z) * mask))
        a = 0.5 + 0.3 * z / (np.std(z) + 1e-12)
    else:
        a = 0.5 + 0.05 * rng.standard_normal(n)
    a = np.clip(a, EPS, 1 - EPS)
    a = a * (n / 2.0) / a.sum()
    return np.clip(a, EPS, 1 - EPS)


def perturb(a, rng, amplitude=0.05, low_freq=True, n_modes=8):
    n = len(a)
    x = np.linspace(-1.0, 1.0, n, endpoint=False) + 1.0 / n
    if low_freq:
        delta = np.zeros(n)
        for k in range(1, n_modes + 1):
            delta += (rng.standard_normal() * np.cos(k * np.pi * x)
                      + rng.standard_normal() * np.sin(k * np.pi * x)) / k
        delta = delta / (np.std(delta) + 1e-12) * amplitude
    else:
        delta = amplitude * rng.standard_normal(n)
    a_new = np.clip(a + delta, EPS, 1 - EPS)
    a_new = a_new * (n / 2.0) / a_new.sum()
    return np.clip(a_new, EPS, 1 - EPS)


def interpolate(a, n_new):
    n_old = len(a)
    x_old = (np.arange(n_old) + 0.5) / n_old
    x_new = (np.arange(n_new) + 0.5) / n_new
    a_new = np.interp(x_new, x_old, a)
    a_new = np.clip(a_new, EPS, 1 - EPS)
    a_new = a_new * (n_new / 2.0) / a_new.sum()
    return np.clip(a_new, EPS, 1 - EPS)


# ---------------------------------------------------------------------------
# Active-set LP-direction polish with exact quadratic line search
# ---------------------------------------------------------------------------

def active_set_polish(a, max_iters=300, tol=1e-12, active_tol=1e-4,
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

        # Gradients of active shifts
        grads = np.zeros((K, n))
        for j, idx in enumerate(active_idx):
            grads[j] = grad_single_shift(a, int(idx - (n - 1)))

        # LP: min t  s.t.  grads @ d - t <= 0;  sum(d)=0;  -1 <= d <= 1
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
        if not res.success or res.x is None:
            break
        d = res.x[:n]; t = res.x[-1]
        if t > -1e-14:
            break

        # Max step from box feasibility
        with np.errstate(divide='ignore', invalid='ignore'):
            pos = d > 1e-15; neg = d < -1e-15
            step_pos = np.min((1.0 - a)[pos] / d[pos]) if pos.any() else np.inf
            step_neg = np.min(a[neg] / (-d[neg])) if neg.any() else np.inf
        max_step = 0.999 * min(step_pos, step_neg)
        if max_step <= 1e-16:
            break

        # Exact quadratic line search: compute lin_k and q_k for all k
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

        def U_of_s(s):
            return float((vals + s * lin + s * s * q).max())

        # Backtracking from max_step
        s = max_step
        best_s, best_u = 0.0, ub
        for _ in range(80):
            u = U_of_s(s)
            if u < best_u - 1e-15:
                best_u, best_s = u, s
                break
            s *= 0.5
            if s < 1e-17:
                break
        # Refine
        if best_s > 0:
            for factor in (1.3, 1.5, 2.0, 0.8, 0.6):
                ss = min(best_s * factor, max_step)
                u = U_of_s(ss)
                if u < best_u - 1e-15:
                    best_u, best_s = u, ss
        if best_s == 0.0:
            break

        a = a + best_s * d
        a = np.clip(a, 0.0, 1.0)
        excess = a.sum() - n / 2.0
        if abs(excess) > 1e-12:
            inner = (a > 1e-6) & (a < 1 - 1e-6)
            if inner.any():
                a[inner] -= excess / inner.sum()
                a = np.clip(a, 0.0, 1.0)

        cur_ub = float(overlap_full(a).max())
        if verbose and it % 20 == 0:
            print(f"  AS iter {it:3d}: K={K:3d}  s={best_s:.3e}  ub={cur_ub:.10f}  d={last_ub - cur_ub:.2e}")
        if last_ub - cur_ub < tol:
            last_ub = cur_ub
            break
        last_ub = cur_ub
    return a


# ---------------------------------------------------------------------------
# Mass-transport hillclimb (SOTA-style)
# Moves: take delta mass from a_i and give to a_j. Preserves sum exactly;
# keeps box constraint exactly when delta <= min(a_i, 1 - a_j).
# delta is a multiple of step (e.g., 2^-32).
# ---------------------------------------------------------------------------

def mass_transport_hillclimb(a, time_budget=300.0, step_init=1e-3,
                             step_min=2.0 ** -32, n_pairs_eval=200,
                             rng_seed=0, verbose=False):
    """Hillclimb with mass-transport moves. At each iteration:
       - compute active shifts (within thr of max)
       - compute gradient direction g (linear LSE-weighted)
       - propose pairs (i, j) prioritizing j with smallest gradient (good place
         to put mass) and i with largest gradient (good place to remove)
       - try each pair with step=delta; accept if max overlap strictly decreases
    Step is halved when stalling, until step_min.
    """
    a = a.copy()
    n = len(a)
    t_start = time.time()
    step = step_init
    cur_ub = float(overlap_full(a).max())
    rng = np.random.default_rng(rng_seed)
    iter_count = 0
    accept_count = 0
    stall_streak = 0

    while time.time() - t_start < time_budget:
        iter_count += 1
        vals = overlap_full(a)
        cur_ub = vals.max()
        thr = cur_ub - 5e-5
        active = np.where(vals >= thr)[0]
        if len(active) == 0:
            break
        # weighted gradient (focus on active shifts)
        w = np.zeros(2 * n - 1)
        ev = np.exp(5000.0 * (vals[active] - cur_ub))
        w[active] = ev / ev.sum()
        g = grad_weighted(a, w)

        # Rank candidate "donor" i (largest g, removing mass helps) and
        # "receiver" j (smallest g, adding mass there helps).
        # Donor must have a_i > 0, receiver must have a_j < 1.
        eligible_donor = np.where(a > step)[0]
        eligible_recv = np.where(a < 1.0 - step)[0]
        if len(eligible_donor) == 0 or len(eligible_recv) == 0:
            step *= 0.5
            if step < step_min:
                break
            continue
        # Top donors and receivers by gradient
        K_top = min(20, len(eligible_donor), len(eligible_recv))
        donor_order = eligible_donor[np.argsort(g[eligible_donor])[::-1][:K_top]]
        recv_order = eligible_recv[np.argsort(g[eligible_recv])[:K_top]]
        # plus a few random eligibles to escape local stalls
        n_rand = min(10, len(eligible_donor), len(eligible_recv))
        donor_extra = rng.choice(eligible_donor, size=n_rand, replace=False)
        recv_extra = rng.choice(eligible_recv, size=n_rand, replace=False)
        donors = np.unique(np.concatenate([donor_order, donor_extra]))
        recvs = np.unique(np.concatenate([recv_order, recv_extra]))

        improved = False
        best_drop = 0.0
        best_pair = None
        best_delta = 0.0

        for i in donors:
            for j in recvs:
                if i == j:
                    continue
                delta = min(step, a[i], 1.0 - a[j])
                if delta <= 0:
                    continue
                # ΔM_k = -δ G[k,i] + δ G[k,j] + (quadratic in δ, only if j = i+k or i = j+k)
                # Compute linear ΔM_k:
                # (grad_k)_i = (2/n) [b[i+k]·1_{0≤i+k≤n-1} - a[i-k]·1_{0≤i-k≤n-1}]
                # Approximate full vals change as δ*(g_j - g_i) per shift (linear).
                # For exactness, also account for quadratic at k = j-i or k = i-j:
                # quad_k = (2/n) * δ^2 * (1[k=j-i] + 1[k=i-j])
                # For small δ this is negligible.
                # Need to compute lin_k = (grad_k)_j - (grad_k)_i for all k.
                # Use per-bin gradient slices.
                lin_k = _per_bin_grad_diff(a, i, j)  # length 2n-1
                d_vals = delta * lin_k
                # Quadratic correction at |k| = |j - i|
                if i != j:
                    k1 = j - i
                    d_vals[k1 + n - 1] += (2.0 / n) * delta * delta
                    k2 = -k1
                    d_vals[k2 + n - 1] += (2.0 / n) * delta * delta
                new_vals_max = float((vals + d_vals).max())
                drop = cur_ub - new_vals_max
                if drop > best_drop + 1e-16:
                    best_drop = drop
                    best_pair = (int(i), int(j))
                    best_delta = float(delta)

        if best_pair is not None and best_drop > 0:
            i, j = best_pair
            a[i] -= best_delta
            a[j] += best_delta
            a[i] = max(0.0, a[i])
            a[j] = min(1.0, a[j])
            accept_count += 1
            stall_streak = 0
            improved = True
            new_ub = float(overlap_full(a).max())
            if verbose and accept_count % 20 == 0:
                print(f"  MT iter {iter_count} accepted {accept_count}, "
                      f"step={step:.2e}, ub={new_ub:.10f}, drop={best_drop:.2e}")
        if not improved:
            stall_streak += 1
            if stall_streak >= 3:
                step *= 0.5
                stall_streak = 0
                if verbose:
                    print(f"  MT step shrunk to {step:.2e}")
                if step < step_min:
                    break
    return a


def _per_bin_grad_diff(a, i, j):
    """Return lin_k = (grad_k)_j - (grad_k)_i for all k in [-(n-1), n-1].
    Length 2n-1, indexed by k+n-1.
    """
    n = len(a)
    b = 1.0 - a
    out = np.zeros(2 * n - 1)
    # (grad_k)_p = (2/n) [b[p+k]·1 - a[p-k]·1] (zero out of range)
    # For each p (=j or i), iterate over k.
    # k+n-1 index.
    # We compute (grad_k)_j - (grad_k)_i across all k.
    # vectorize: for p=j and p=i, the support set of k for which b[p+k] is valid
    # is k in [-p, n-1-p]; similarly for a[p-k]: k in [p-(n-1), p].
    # Build full arrays.
    def bin_grad(p):
        g = np.zeros(2 * n - 1)
        # b-term: g[k+n-1] += b[p+k] for k in [-p, n-1-p]
        # k = m - p, m in [0, n-1], so k+n-1 = m - p + n - 1
        idx = np.arange(n) - p + (n - 1)
        g[idx] += b
        # a-term: g[k+n-1] -= a[p-k] for k in [p-(n-1), p]
        # k = p - m, m in [0, n-1], so k+n-1 = p - m + n - 1
        idx2 = p - np.arange(n) + (n - 1)
        g[idx2] -= a
        return g * (2.0 / n)
    return bin_grad(j) - bin_grad(i)


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

def solve(verbose=False, time_budget_sec=900, target_n=600):
    here = os.path.dirname(os.path.abspath(__file__))
    t_start = time.time()

    init_kinds = ["bell", "cosine", "double_cos", "triple_cos",
                  "tent", "raised_cos", "two_block",
                  "rand_smooth", "bell_narrow", "bell_wide"]
    betas_main = (50.0, 200.0, 800.0, 3000.0, 10000.0)
    betas_polish = (3000.0, 15000.0, 50000.0)
    iters_main = 300

    # ---- Phase A: multi-init at n=target_n ----
    n = target_n
    best_a, best_ub = None, np.inf
    for ki, kind in enumerate(init_kinds):
        if time.time() - t_start > time_budget_sec * 0.25:
            break
        a0 = make_init(n, kind, np.random.default_rng(11 + ki))
        a_opt, ub = optimize_from_a(a0, betas_main, iters_main)
        if verbose:
            print(f"[init n={n} {kind:<14}] ub={ub:.6f}")
        if ub < best_ub:
            best_a, best_ub = a_opt, ub
    if verbose:
        print(f"[after init] best ub = {best_ub:.6f}  t={time.time()-t_start:.0f}s")

    # ---- Phase B: basin hopping at n=target_n ----
    rng = np.random.default_rng(101)
    hop = 0
    while time.time() - t_start < time_budget_sec * 0.50:
        hop += 1
        amp = float(rng.choice([0.02, 0.05, 0.10, 0.20]))
        lf = bool(rng.integers(0, 2))
        modes = int(rng.integers(3, 20))
        a_try = perturb(best_a, rng, amplitude=amp, low_freq=lf, n_modes=modes)
        a_opt, ub = optimize_from_a(a_try, (200.0, 1000.0, 5000.0, 30000.0), 200)
        if ub < best_ub - 1e-7:
            if verbose:
                print(f"[hop #{hop} amp={amp} lf={lf} modes={modes}] ub={ub:.6f} (NEW)")
            best_a, best_ub = a_opt, ub

    # ---- Phase C: LSE polish at high beta ----
    a_opt, ub = optimize_from_a(best_a, betas_polish, 600)
    if ub < best_ub:
        if verbose: print(f"[LSE polish] ub={ub:.6f} (NEW)")
        best_a, best_ub = a_opt, ub

    # ---- Phase D: active-set LP polish ----
    polished = active_set_polish(best_a, max_iters=400, tol=1e-13,
                                 active_tol=1e-4, verbose=verbose)
    ub_p = float(overlap_full(polished).max())
    if ub_p < best_ub:
        if verbose: print(f"[active-set polish] ub={ub_p:.10f} (NEW)")
        best_a, best_ub = polished, ub_p

    # ---- Phase E: mass-transport hillclimb (SOTA-style) ----
    remaining = max(60.0, time_budget_sec - (time.time() - t_start))
    if verbose:
        print(f"[mass transport] starting; budget {remaining:.0f}s, "
              f"current ub = {best_ub:.10f}")
    a_mt = mass_transport_hillclimb(best_a, time_budget=remaining,
                                    step_init=1e-3, step_min=2.0 ** -32,
                                    verbose=verbose)
    ub_mt = float(overlap_full(a_mt).max())
    if ub_mt < best_ub:
        if verbose: print(f"[mass transport] final ub={ub_mt:.10f} (NEW)")
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
    target_n = int(os.environ.get("TARGET_N", "600"))
    print(json.dumps(solve(verbose=True, time_budget_sec=budget,
                           target_n=target_n), indent=2))
