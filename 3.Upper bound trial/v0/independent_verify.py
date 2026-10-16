"""
Adversarial *independent* verifier for the proposal-agent sequence.

Goals (per verifier.md):
  - Recompute U(a) = max_k M_k(a) with implementations NOT shared
    with the proposer's overlap_full / overlap_fft routines.
  - Detect indexing bugs, periodicity mistakes, normalization errors,
    hidden infeasibility, discretization artifacts, missed active shifts.
  - Report certified value, active shifts, near-active shifts, feasibility
    residuals, refinement / Fourier diagnostics, reproducibility metadata.

Independent computations used:
  (B) pure-Python double loop over (k, i)           -- ground truth
  (C) np.einsum-based matrix construction           -- second path
  (D) shift-and-multiply (slicing) per k            -- third path
  (E) zero-pad CONVOLUTION via np.convolve, then
      reflect-and-reindex                            -- fourth path
  (F) 4x piecewise-constant refinement of a, FFT eval
      to check that finer-grid sup is no smaller

All four (B,C,D,E) must agree to ~1e-12 and must reproduce the claim.
"""
from __future__ import annotations
import json, os, sys, hashlib, platform, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))


# ---------- four independent recomputations of M_k --------------------------

def overlap_loop(a: np.ndarray) -> np.ndarray:
    """B: pure Python double loop. Slow but unambiguous ground truth."""
    n = a.size
    out = np.zeros(2 * n - 1, dtype=np.float64)
    for k in range(-(n - 1), n):
        i_lo = max(0, -k)
        i_hi = min(n, n - k)
        acc = 0.0
        for i in range(i_lo, i_hi):
            acc += a[i] * (1.0 - a[i + k])
        out[k + (n - 1)] = (2.0 / n) * acc
    return out


def overlap_einsum(a: np.ndarray) -> np.ndarray:
    """C: build outer-product matrix and pick out diagonals via einsum."""
    n = a.size
    b = 1.0 - a
    P = np.einsum("i,j->ij", a, b)            # P[i,j] = a_i * (1 - a_j)
    out = np.empty(2 * n - 1)
    # diagonal offset k: i,j=i+k; sum over i in [max(0,-k), min(n,n-k))
    for k in range(-(n - 1), n):
        out[k + (n - 1)] = P.diagonal(offset=k).sum()
    return out * (2.0 / n)


def overlap_slice(a: np.ndarray) -> np.ndarray:
    """D: shift-and-multiply per k using numpy slicing only."""
    n = a.size
    b = 1.0 - a
    out = np.empty(2 * n - 1)
    for k in range(-(n - 1), n):
        if k >= 0:
            out[k + (n - 1)] = float(np.dot(a[: n - k], b[k:]))
        else:
            kk = -k
            out[k + (n - 1)] = float(np.dot(a[kk:], b[: n - kk]))
    return out * (2.0 / n)


def overlap_convolve(a: np.ndarray) -> np.ndarray:
    """E: get M_k via np.convolve(a, reverse(b)).

    Convolution:  c[m] = sum_j a[j] * br[m - j]   where br = b[::-1].
    With m = i + (n-1-(i+k)) = n-1-k, we get c[n-1-k] = sum_i a[i] b[i+k].
    So  M_k = (2/n) * c[n-1-k].
    """
    n = a.size
    b = 1.0 - a
    c = np.convolve(a, b[::-1], mode="full")        # length 2n-1
    out = np.empty(2 * n - 1)
    for k in range(-(n - 1), n):
        out[k + (n - 1)] = c[n - 1 - k]
    return out * (2.0 / n)


# ---------- piecewise-constant refinement ------------------------------------

def refine_pwc(a: np.ndarray, factor: int) -> np.ndarray:
    """Each bin -> `factor` identical sub-bins; encodes the SAME continuous
    step function, but on a finer mesh, so refined sup is honest."""
    return np.repeat(a, factor)


def overlap_fft_indep(a: np.ndarray) -> np.ndarray:
    """Indep. FFT recompute used only for the refined grid (where loops
    are too slow). Different padding strategy than the proposer's."""
    n = a.size
    b = 1.0 - a
    # zero-pad to exactly 2n+1 (not next power of two) -> different code path
    Lpad = 2 * n + 1
    A = np.fft.fft(np.concatenate([a, np.zeros(Lpad - n)]))
    B = np.fft.fft(np.concatenate([b, np.zeros(Lpad - n)]))
    full = np.real(np.fft.ifft(np.conjugate(A) * B))
    out = np.empty(2 * n - 1)
    for k in range(-(n - 1), n):
        out[k + (n - 1)] = full[k % Lpad]
    return out * (2.0 / n)


# ---------- main -------------------------------------------------------------

def verify(seq_path: str, max_loop_n: int = 1500) -> dict:
    with open(seq_path, "rb") as f:
        raw = f.read()
    digest = hashlib.sha256(raw).hexdigest()
    data = json.loads(raw)
    a = np.asarray(data["a"], dtype=np.float64)
    n = a.size
    claimed_ub = float(data["upper_bound"])
    claimed_k = int(data["argmax_shift"])
    claimed_n = int(data.get("n", n))

    # ---- feasibility --------------------------------------------------------
    box_lo_res = float(min(0.0, a.min()))
    box_hi_res = float(max(0.0, a.max() - 1.0))
    mass_err   = float(abs(a.sum() - n / 2.0))
    feas_ok    = (box_lo_res >= -1e-9 and box_hi_res <= 1e-9 and mass_err <= 1e-8)

    # ---- four independent recomputations -----------------------------------
    t0 = time.time()
    if n <= max_loop_n:
        vals_loop = overlap_loop(a)
    else:
        vals_loop = None
    t_loop = time.time() - t0

    t0 = time.time();  vals_es = overlap_einsum(a)   ;  t_es   = time.time() - t0
    t0 = time.time();  vals_sl = overlap_slice(a)    ;  t_sl   = time.time() - t0
    t0 = time.time();  vals_cv = overlap_convolve(a) ;  t_cv   = time.time() - t0

    if vals_loop is not None:
        d_loop_es = float(np.max(np.abs(vals_loop - vals_es)))
        d_loop_sl = float(np.max(np.abs(vals_loop - vals_sl)))
        d_loop_cv = float(np.max(np.abs(vals_loop - vals_cv)))
        truth = vals_loop
        truth_name = "loop"
    else:
        d_loop_es = d_loop_sl = d_loop_cv = float("nan")
        truth = vals_es
        truth_name = "einsum"
    d_es_sl = float(np.max(np.abs(vals_es - vals_sl)))
    d_es_cv = float(np.max(np.abs(vals_es - vals_cv)))
    d_sl_cv = float(np.max(np.abs(vals_sl - vals_cv)))

    ub_truth   = float(truth.max())
    k_truth    = int(np.argmax(truth) - (n - 1))
    ub_diff    = float(abs(ub_truth - claimed_ub))
    k_match    = (k_truth == claimed_k)

    # ---- near-active shifts -------------------------------------------------
    order = np.argsort(truth)[::-1]
    near = []
    for idx in order[:15]:
        near.append({
            "k": int(idx - (n - 1)),
            "M_k": float(truth[idx]),
            "gap_to_max": float(ub_truth - truth[idx]),
        })

    # ---- check missed active shifts at boundary |k| ≈ n --------------------
    boundary_max = float(max(truth[0], truth[-1]))

    # ---- refinement to finer grid ------------------------------------------
    factor = 4
    a_fine = refine_pwc(a, factor)
    vals_fine = overlap_fft_indep(a_fine)
    ub_fine   = float(vals_fine.max())
    k_fine    = int(np.argmax(vals_fine) - (a_fine.size - 1))
    refine_ok = (ub_fine <= ub_truth + 1e-10)

    # ---- Fourier energy concentration -------------------------------------
    A_hat = np.fft.rfft(a - a.mean())
    energy = np.abs(A_hat) ** 2
    energy = energy / max(energy.sum(), 1e-300)
    top_freqs = np.argsort(energy)[::-1][:5]
    fourier = [
        {"freq_bin": int(f), "rel_energy": float(energy[f])} for f in top_freqs
    ]

    # ---- value / impl tolerances ------------------------------------------
    impl_ok = (max(d_es_sl, d_es_cv, d_sl_cv,
                   d_loop_es if vals_loop is not None else 0.0,
                   d_loop_sl if vals_loop is not None else 0.0,
                   d_loop_cv if vals_loop is not None else 0.0) <= 1e-11)
    value_ok = (ub_diff <= 1e-10) and k_match
    verdict  = "VALID" if (feas_ok and value_ok and impl_ok and refine_ok) else "INVALID"

    return {
        "verdict": verdict,
        "certified_upper_bound": ub_truth,
        "truth_implementation": truth_name,
        "claimed_upper_bound": claimed_ub,
        "ub_discrepancy": ub_diff,
        "argmax_k_truth": k_truth,
        "argmax_k_claimed": claimed_k,
        "argmax_shift_matches": k_match,
        "implementation_cross_check": {
            "loop_vs_einsum":  d_loop_es,
            "loop_vs_slice":   d_loop_sl,
            "loop_vs_convol":  d_loop_cv,
            "einsum_vs_slice": d_es_sl,
            "einsum_vs_convol": d_es_cv,
            "slice_vs_convol": d_sl_cv,
            "ok": impl_ok,
        },
        "timings_sec": {
            "loop": t_loop if vals_loop is not None else None,
            "einsum": t_es, "slice": t_sl, "convolve": t_cv,
        },
        "feasibility": {
            "box_lower_residual": box_lo_res,
            "box_upper_residual": box_hi_res,
            "mass_error": mass_err,
            "n_matches_file": (n == claimed_n),
            "ok": feas_ok,
        },
        "boundary_shift_max_M": boundary_max,
        "refinement": {
            "factor": factor,
            "ub_on_fine_grid": ub_fine,
            "argmax_k_on_fine_grid": k_fine,
            "monotone_ok": refine_ok,
        },
        "near_active_shifts": near,
        "fourier_top_modes": fourier,
        "metadata": {
            "n": n,
            "sequence_sha256": digest,
            "python": sys.version.split()[0],
            "numpy": np.__version__,
            "platform": platform.platform(),
            "verifier_file": __file__,
        },
    }


if __name__ == "__main__":
    seq = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "sequence.json")
    report = verify(seq)
    out_path = os.path.join(HERE, "independent_verifier_report.json")
    with open(out_path, "w") as f:
        json.dump(report, f, indent=2)
    print(json.dumps(report, indent=2))
