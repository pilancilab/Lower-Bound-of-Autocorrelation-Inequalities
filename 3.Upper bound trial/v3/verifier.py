"""
Independent verifier for the proposed sequence (per verifier.md).

Recomputes M_k via:
  - O(n^2) naive double loop (reference)
  - FFT-based linear cross-correlation (fast, cross-check)

Checks feasibility (box, mass), reports active argmax shift,
near-active shifts (within 1e-4 of max), and FFT-vs-naive
maximum absolute discrepancy.
"""

import json
import os
import sys
import hashlib
import numpy as np


def overlap_naive(a):
    n = len(a)
    M = np.zeros(2 * n - 1)
    for d in range(-(n - 1), n):
        i_lo = max(0, -d)
        i_hi = min(n, n - d)
        s = 0.0
        for i in range(i_lo, i_hi):
            s += a[i] * (1.0 - a[i + d])
        M[d + n - 1] = (2.0 / n) * s
    return M


def overlap_fft(a):
    n = len(a)
    b = 1.0 - a
    L = 1
    while L < 2 * n:
        L *= 2
    A = np.fft.rfft(a, n=L)
    B = np.fft.rfft(b, n=L)
    full = np.fft.irfft(np.conjugate(A) * B, n=L)
    M = np.zeros(2 * n - 1)
    for d in range(-(n - 1), n):
        M[d + n - 1] = (2.0 / n) * full[d % L]
    return M


def verify(seq_path, naive_n_limit=2000):
    with open(seq_path) as f:
        data = json.load(f)
    a = np.asarray(data["a"], dtype=float)
    n = len(a)
    claimed_n = int(data.get("n", n))
    claimed_ub = float(data.get("upper_bound", float("nan")))
    claimed_arg = data.get("argmax_shift", None)

    box_lo = float(a.min())
    box_hi = float(a.max())
    mass = float(a.sum())
    mass_target = n / 2.0
    mass_err = abs(mass - mass_target)

    if box_lo < -1e-9 or box_hi > 1 + 1e-9:
        feas = "BOX_VIOLATION"
    elif mass_err > 1e-6 * max(1.0, mass_target):
        feas = "MASS_VIOLATION"
    else:
        feas = "ok"

    use_naive = n <= naive_n_limit
    if use_naive:
        M_naive = overlap_naive(a)
    M_fft = overlap_fft(a)
    cross_err = (float(np.max(np.abs(M_naive - M_fft))) if use_naive
                 else float("nan"))

    ub = float(M_fft.max())
    arg = int(np.argmax(M_fft) - (n - 1))

    tol = 1e-4
    near = [(int(i - (n - 1)), float(M_fft[i]))
            for i in range(2 * n - 1) if M_fft[i] >= ub - tol]
    near.sort(key=lambda t: -t[1])
    sha = hashlib.sha256(a.tobytes()).hexdigest()[:16]

    print("=" * 64)
    print(f"Sequence:           {seq_path}")
    print(f"n (file):           {claimed_n}  (len(a)={n})")
    print(f"sha256[:16](a):     {sha}")
    print(f"Box range:          [{box_lo:.3e}, {box_hi:.6f}]")
    print(f"Sum a:              {mass:.10f}  (target {mass_target:.10f}, err {mass_err:.2e})")
    print(f"Feasibility:        {feas}")
    print("-" * 64)
    print(f"Certified U(a) (FFT):     {ub:.12f}  at shift d = {arg}")
    if use_naive:
        print(f"Naive recomputation:      {float(M_naive.max()):.12f}  at d = {int(np.argmax(M_naive)-(n-1))}")
        print(f"FFT vs naive max abs err: {cross_err:.2e}")
    print(f"Claimed upper_bound:      {claimed_ub:.12f}")
    print(f"Claimed argmax_shift:     {claimed_arg}")
    diff = ub - claimed_ub
    print(f"Verifier - claimed:       {diff:+.2e}")
    print("-" * 64)
    print(f"Near-active shifts within {tol}:")
    for d, v in near[:20]:
        print(f"    d = {d:+5d}    M_d = {v:.10f}")
    print("=" * 64)

    valid = (feas == "ok"
             and (np.isnan(cross_err) or cross_err < 1e-9)
             and abs(diff) < 1e-6)
    print(f"VERDICT: {'VALID' if valid else 'INVALID'}")
    return {
        "valid": bool(valid),
        "certified_ub": ub,
        "argmax_shift": arg,
        "feasibility": feas,
        "mass_error": mass_err,
        "box_lo": box_lo, "box_hi": box_hi,
        "cross_check_err": cross_err,
        "claimed_ub": claimed_ub,
        "near_active": near[:20],
        "sha": sha,
    }


if __name__ == "__main__":
    p = sys.argv[1] if len(sys.argv) > 1 else os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "sequence.json"
    )
    verify(p)
