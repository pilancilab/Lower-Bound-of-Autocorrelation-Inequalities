"""
v9 — CMA-ES on asymmetric K-level step functions, then full polish.

We parameterize f as K constant levels on K subintervals of [-1, 1] with
K-1 internal breakpoints. Asymmetric (no enforced symmetry). Free
parameters: K heights in [0, 1] and K-1 breakpoints; mass is enforced by
rescaling the heights to sum to n/2 in discrete form.

For each K in {4, 5, 6, 7, 8, 10}, run CMA-ES with multiple restarts.
Take the best parametric solution, convert to n=600, polish with the full
pipeline (LSE → AS → MT).
"""

import json, os, time
import numpy as np
from scipy.optimize import minimize, brentq, linprog
import cma

EPS = 1e-9


# ---- core utilities -------------------------------------------------------

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


# ---- asymmetric K-level encoding ------------------------------------------

def theta_to_a(theta, K, n):
    """theta is a 2K-1 vector: theta[:K] -> heights in (0,1) via sigmoid;
    theta[K:] -> K-1 internal breakpoint params -> breakpoints in (-1, 1).
    """
    raw_h = theta[:K]
    raw_b = theta[K:]
    # Heights via sigmoid
    h = 1.0 / (1.0 + np.exp(-raw_h))
    # Breakpoints: use cumulative-softmax to ensure ordering in (-1, 1).
    # K-1 breakpoints partition [-1, 1] into K intervals. Use widths.
    # widths = softmax(raw_b ++ a learned dummy), then breakpoints = -1 + 2 * cumsum(widths[:K-1])
    # For 2K-1 params total, raw_b has K-1 entries. We'll create K widths
    # by interpreting raw_b as logits for K-1 widths and the K-th is 1 - sum.
    # Actually simplest: include a fixed dummy 0 as the K-th width logit.
    logits = np.concatenate([raw_b, [0.0]])
    e = np.exp(logits - logits.max())
    widths = e / e.sum()
    # Cumulative on [0, 1], then map to [-1, 1]
    cum = np.cumsum(widths)[:-1]  # K-1 breakpoints in (0, 1)
    breaks = 2 * cum - 1  # in (-1, 1)
    # Build a on n bins
    x_centers = np.linspace(-1.0, 1.0, n, endpoint=False) + 1.0 / n
    # Determine which interval each x_center falls in
    bidx = np.searchsorted(breaks, x_centers, side='right')  # in [0, K-1]
    a = h[bidx]
    return normalize(a)


def klevel_obj(theta, K, n_eval):
    a = theta_to_a(theta, K, n_eval)
    return float(overlap_full(a).max())


def klevel_cma(K, n_eval, time_budget, sigma=1.0, seed=0, verbose=False):
    """Run CMA-ES on the K-level parametric family with multiple restarts."""
    n_param = 2 * K - 1
    t_start = time.time()
    best_theta, best_ub = None, np.inf
    rng = np.random.default_rng(seed)
    restart = 0
    while time.time() - t_start < time_budget:
        restart += 1
        # Random initial theta in [-3, 3]
        theta0 = rng.uniform(-3.0, 3.0, size=n_param)
        try:
            es = cma.CMAEvolutionStrategy(
                theta0, sigma,
                {'verbose': -9, 'maxiter': 80, 'tolfun': 1e-9,
                 'tolx': 1e-7, 'popsize': max(8, 4 + int(3 * np.log(n_param)))}
            )
            while not es.stop() and time.time() - t_start < time_budget:
                solutions = es.ask()
                fitnesses = [klevel_obj(s, K, n_eval) for s in solutions]
                es.tell(solutions, fitnesses)
            x_best = es.result.xbest
            f_best = es.result.fbest
            if f_best is not None and f_best < best_ub:
                best_ub = float(f_best)
                best_theta = np.asarray(x_best)
                if verbose:
                    print(f"   K={K} restart {restart}: ub={best_ub:.10f}")
        except Exception as e:
            if verbose: print(f"   K={K} restart {restart} failed: {e}")
            continue
    return best_theta, best_ub


# ---- polish ----------------------------------------------------------------

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
                          step_min=2.0 ** -38, verbose=False, log_every=200):
    a = a.copy(); n = len(a); t_start = time.time()
    step = step_init; accept = 0; stall = 0; itr = 0
    rng = np.random.default_rng(777)
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
        K_top = min(30, len(elig_d), len(elig_r))
        donor_top = elig_d[np.argsort(g[elig_d])[::-1][:K_top]]
        recv_top = elig_r[np.argsort(g[elig_r])[:K_top]]
        donor_rand = rng.choice(elig_d, size=min(15, len(elig_d)), replace=False)
        recv_rand = rng.choice(elig_r, size=min(15, len(elig_r)), replace=False)
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


# ---- driver ---------------------------------------------------------------

def solve(verbose=False, time_budget_sec=900):
    here = os.path.dirname(os.path.abspath(__file__))
    t_start = time.time()
    n_main = 600

    # Warm from v8
    v8_path = os.path.join(here, "..", "v8", "sequence.json")
    with open(v8_path) as f:
        v8_data = json.load(f)
    a_v8 = np.asarray(v8_data["a"], dtype=float)
    ub_v8 = float(v8_data.get("upper_bound", float(overlap_full(a_v8).max())))
    if verbose:
        print(f"[warm v8] ub={ub_v8:.17f}")
    best_a, best_ub = a_v8, ub_v8

    # ---- Phase 1: CMA-ES on asymmetric K-level for K in {4..10} ----
    K_budget = max(20.0, time_budget_sec * 0.6 / 5)  # split across K values
    n_eval = 300  # coarse eval for speed
    klevel_results = []
    for K in (4, 5, 6, 8, 10):
        if time.time() - t_start > time_budget_sec * 0.65: break
        if verbose: print(f"\n[CMA-ES K={K}, n_eval={n_eval}, budget {K_budget:.0f}s]")
        theta, ub = klevel_cma(K, n_eval, K_budget, sigma=1.5,
                               seed=42 + K, verbose=verbose)
        if theta is not None:
            if verbose: print(f"   K={K} CMA-ES best @ n=300: {ub:.10f}")
            # Convert to n=600 and polish
            a_init = theta_to_a(theta, K, n_main)
            a_opt, ub600 = optimize_from_a(a_init, (200, 1500, 10000, 50000), 250)
            if verbose: print(f"   -> n=600 LBFGS: {ub600:.10f}")
            if ub600 < best_ub - 1e-9:
                if verbose: print(f"   *** K={K} NEW BEST: {ub600:.17f}")
                best_a, best_ub = a_opt, ub600

    # ---- Phase 2: final polish ----
    if verbose: print(f"\n[Final polish] current {best_ub:.17f}")
    a_opt, ub = optimize_from_a(best_a, (5000, 25000, 100000), 500)
    if ub < best_ub:
        if verbose: print(f"[LSE] {ub:.17f} (NEW)")
        best_a, best_ub = a_opt, ub
    pol = active_set_polish(best_a, max_iters=300, tol=1e-14)
    ub_p = float(overlap_full(pol).max())
    if ub_p < best_ub:
        if verbose: print(f"[AS] {ub_p:.17f} (NEW)")
        best_a, best_ub = pol, ub_p
    remaining = max(45.0, time_budget_sec - (time.time() - t_start))
    if verbose: print(f"[MT] budget {remaining:.0f}s, cur {best_ub:.17f}")
    a_mt = mass_transport_polish(best_a, time_budget=remaining,
                                 step_init=2.0 ** -16, step_min=2.0 ** -38,
                                 verbose=verbose)
    ub_mt = float(overlap_full(a_mt).max())
    if ub_mt < best_ub:
        if verbose: print(f"[MT] {ub_mt:.17f} (NEW)")
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
