"""
v15 — Exact integer-arithmetic pair-move search.

v13/14 plateaued because float64 noise (~1e-13) hides true descent at this
basin's local minimum. Exact arithmetic reveals that only shift d = +220
is the true maximum; many others appear tied in float64 but are strictly
below at machine precision.

Strategy:
  - Quantize a to integers p on a 2^-Q dyadic grid (Q = 48).
  - Maintain cached S_k = sum_i p_i (Q_int - p_{i+k}) for all k as Python ints.
  - For a pair move (i, j, δ=1) [smallest possible], compute ΔS_k for all
    shifts via a closed-form O(n) update of S, then evaluate the new max
    in exact arithmetic.
  - Accept any move that strictly decreases the integer max.
  - Iterate until no improving (i, j) found among top candidates ranked by
    a float64 surrogate.
"""

import json, os, time, sys
import numpy as np

EPS = 1e-9
QBITS = 48
Q_INT = 1 << QBITS  # the dyadic denominator

HERE = os.path.dirname(os.path.abspath(__file__))


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


def quantize_and_repair(a):
    """Quantize a to integer multiples of 1/Q_INT, ensure sum = n * Q_INT // 2 exactly."""
    n = len(a)
    p = [int(round(x * Q_INT)) for x in a]
    p = [max(0, min(Q_INT, x)) for x in p]
    target = n * Q_INT // 2
    diff = sum(p) - target
    # Repair: bump entries by ±1
    if diff != 0:
        order = list(range(n))
        # prefer interior entries
        order.sort(key=lambda i: -(min(p[i], Q_INT - p[i])))
        sign = -1 if diff > 0 else 1
        for idx in order:
            if diff == 0: break
            new = p[idx] + sign
            if 0 <= new <= Q_INT:
                p[idx] = new
                diff += sign
    assert sum(p) == target, f"repair failed: {sum(p)} != {target}"
    return p


def compute_all_Sk(p):
    """Returns list S of length 2n-1, S[d+n-1] = sum_i p_i (Q_INT - p_{i+d}).

    We compute by direct loop. O(n^2). For n=600 takes ~1s in pure Python.
    """
    n = len(p)
    S = [0] * (2 * n - 1)
    for d in range(-(n - 1), n):
        s = 0
        if d >= 0:
            for i in range(n - d):
                s += p[i] * (Q_INT - p[i + d])
        else:
            kk = -d
            for i in range(kk, n):
                s += p[i] * (Q_INT - p[i - kk])
        S[d + n - 1] = s
    return S


def update_Sk_after_pair_move(S, p, i, j, n):
    """Compute new_S after move p[i] -= 1, p[j] += 1, given current S and p.
    Returns new_S as a NEW list (does not modify input).
    """
    new_S = S.copy()
    Q = Q_INT
    for d in range(-(n - 1), n):
        d_idx = d + n - 1
        delta = 0
        # l = i contribution
        if 0 <= i + d < n:
            # old: p[i] * (Q - p[i+d])
            # new: (p[i]-1) * (Q - p[i+d] - (-1 if i+d == j else 0)) if i+d == j then p_new[i+d] = p[j]+1
            pik = p[i + d]
            if i + d == j: pik_new = p[j] + 1
            else: pik_new = pik
            delta += (p[i] - 1) * (Q - pik_new) - p[i] * (Q - pik)
        # l = j contribution
        if 0 <= j + d < n:
            pjk = p[j + d]
            if j + d == i: pjk_new = p[i] - 1
            else: pjk_new = pjk
            delta += (p[j] + 1) * (Q - pjk_new) - p[j] * (Q - pjk)
        # l = i - d contribution (where l+d = i)
        l = i - d
        if 0 <= l < n and l != i and l != j:
            # p_l unchanged, (Q - p[i]) -> (Q - (p[i]-1)) = (Q - p[i]) + 1
            delta += p[l] * 1
        # l = j - d contribution (where l+d = j)
        l = j - d
        if 0 <= l < n and l != i and l != j:
            # p_l unchanged, (Q - p[j]) -> (Q - (p[j]+1)) = (Q - p[j]) - 1
            delta -= p[l]
        new_S[d_idx] += delta
    return new_S


def search_pair_moves(a, S, p, n_top=80, time_budget=120.0, verbose=False):
    """Hill-climb with exact integer arithmetic. For δ_int=1 only.

    Returns (final p, final S, accepted count).
    """
    n = len(p)
    t_start = time.time()
    accept = 0
    iters = 0
    rng = np.random.default_rng(42)
    while time.time() - t_start < time_budget:
        iters += 1
        cur_max = max(S)
        cur_argmax = S.index(cur_max)
        cur_d = cur_argmax - (n - 1)
        # Gradient w.r.t. argmax shift in float
        vals = overlap_full(np.asarray(a, dtype=float))
        # Pick donors and receivers via gradient of weighted active shifts
        ub_f = vals.max()
        active = np.where(vals >= ub_f - 1e-7)[0]
        if len(active) == 0: active = np.array([cur_argmax])
        w = np.zeros(2 * n - 1)
        w[active] = np.exp(5000.0 * (vals[active] - ub_f))
        w /= w.sum()
        g = grad_weighted(np.asarray(a, dtype=float), w)
        # Eligible donors (p[i] > 0) and receivers (p[j] < Q_INT)
        elig_d = [i for i in range(n) if p[i] > 0]
        elig_r = [j for j in range(n) if p[j] < Q_INT]
        if not elig_d or not elig_r: break
        # Rank
        donors = sorted(elig_d, key=lambda i: -g[i])[:n_top]
        recvs = sorted(elig_r, key=lambda j: g[j])[:n_top]
        # Plus random
        donors_rand = list(rng.choice(elig_d, size=min(40, len(elig_d)), replace=False))
        recvs_rand = list(rng.choice(elig_r, size=min(40, len(elig_r)), replace=False))
        candidates = set()
        for di in donors + donors_rand:
            for re in recvs + recvs_rand:
                if di != re: candidates.add((di, re))

        best_max = cur_max
        best_pair = None
        best_S = None
        for (di, re) in candidates:
            new_S = update_Sk_after_pair_move(S, p, di, re, n)
            new_max = max(new_S)
            if new_max < best_max:
                best_max = new_max
                best_pair = (di, re)
                best_S = new_S
        if best_pair is None:
            if verbose: print(f"  iter {iters}: no improving move among {len(candidates)} candidates")
            break
        di, re = best_pair
        p[di] -= 1; p[re] += 1
        a[di] = p[di] / Q_INT
        a[re] = p[re] / Q_INT
        S = best_S
        accept += 1
        if verbose and accept % 10 == 0:
            cur_max_real = best_max * 2.0 / (n * Q_INT * Q_INT)
            print(f"  iter {iters} accepts {accept} move ({di},{re}) "
                  f"new_max_int={best_max} ub={cur_max_real:.20f}")
    return p, a, S, accept


def solve(verbose=False, time_budget_sec=600):
    here = os.path.dirname(os.path.abspath(__file__))
    t_start = time.time()

    # Load v13
    v13_path = os.path.join(here, "..", "v13", "sequence.json")
    with open(v13_path) as f:
        d = json.load(f)
    a = list(np.asarray(d["a"], dtype=float))
    n = len(a)
    ub_float = float(overlap_full(np.asarray(a)).max())
    if verbose:
        print(f"[warm v13] float U = {ub_float:.17f}, n = {n}")

    # Quantize to 2^-48 grid
    p = quantize_and_repair(a)
    a = [pi / Q_INT for pi in p]
    if verbose: print(f"[quantized to 2^-{QBITS}] sum(p)={sum(p)}, target={n*Q_INT//2}")

    # Compute initial S
    t0 = time.time()
    S = compute_all_Sk(p)
    if verbose: print(f"[initial S] computed in {time.time()-t0:.1f}s, max_int={max(S)}, "
                      f"max as float = {max(S)*2.0/(n*Q_INT*Q_INT):.20f}")

    # Hill climb with exact integer arithmetic
    if verbose: print(f"\n[exact-int pair search] budget {time_budget_sec*0.9:.0f}s")
    p, a, S, accept = search_pair_moves(a, S, p,
                                         n_top=80,
                                         time_budget=time_budget_sec * 0.9,
                                         verbose=verbose)
    if verbose: print(f"\n[done] accepted {accept} moves")

    best_ub_int = max(S)
    best_ub = best_ub_int * 2.0 / (n * Q_INT * Q_INT)
    a_arr = np.asarray(a, dtype=np.float64)
    ub_recomputed = float(overlap_full(a_arr).max())
    if verbose: print(f"[final] integer-exact max as float: {best_ub:.20f}")
    if verbose: print(f"[final] float64 FFT recomputation:  {ub_recomputed:.17f}")

    vals = overlap_full(a_arr)
    arg_d = int(np.argmax(vals) - (n - 1))
    seq_path = os.path.join(here, "sequence.json")
    with open(seq_path, "w") as f:
        json.dump({"n": n, "a": a_arr.tolist(),
                   "upper_bound": float(ub_recomputed),
                   "argmax_shift": arg_d,
                   "upper_bound_exact_str": str(best_ub_int) + " / " + str(n * Q_INT * Q_INT // 2)}, f)
    return {
        "upper_bound": float(ub_recomputed),
        "status": "ok",
        "sequence_path": seq_path,
        "n": int(n),
        "mass_error": float(abs(a_arr.sum() - n / 2.0)),
        "box_violation": float(max(0.0, a_arr.max() - 1.0, -a_arr.min())),
        "argmax_shift": arg_d,
    }


if __name__ == "__main__":
    budget = int(os.environ.get("BUDGET_SEC", "600"))
    print(json.dumps(solve(verbose=True, time_budget_sec=budget), indent=2))
