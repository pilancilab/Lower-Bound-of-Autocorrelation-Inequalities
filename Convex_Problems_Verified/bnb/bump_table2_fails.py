"""Bump N/T on the 5 Table 2 rows that failed at the paper's resolution.

Failures from verify_table2.py (Omega < 0.38):
  A6  N=10k T=4k   h(0,0.08) p(0,1)        q(-1,    -0.05)    Omega=0.37716
  A7  N=10k T=4k   h(0,0.08) p(0,1)        q(-0.05, -0.025)   Omega=0.37451
  A9  N=10k T=4k   h(0,0.08) p(0,1)        q( 0.025, 0.05)    Omega=0.37594
  B1  N=20k T=5k   h(0,0.06) p(0.33,0.45)  q(-0.025,-0.020)   Omega=0.37666
  B2  N=20k T=5k   h(0,0.06) p(0.33,0.45)  q( 0.020, 0.025)   Omega=0.37791

Batch 1: A6, A7, A9 at N=20000, T=5000  (one tier up).
Batch 2: ALL 5  at N=50000, T=10000     (max bump).
"""
from __future__ import annotations

import os
os.environ.setdefault("OMP_NUM_THREADS", "8")
os.environ.setdefault("MKL_NUM_THREADS", "8")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "8")

import json
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT))

import cvxpy as cp
from EMOP import WhiteConvexProblem

THRESHOLD = 0.38

# (label, h1, h2, p1, p2, q1, q2)
FAILS = [
    ("A6", 0.00, 0.08, 0.00, 1.00, -1.000,  -0.050),
    ("A7", 0.00, 0.08, 0.00, 1.00, -0.050,  -0.025),
    ("A9", 0.00, 0.08, 0.00, 1.00,  0.025,   0.050),
    ("B1", 0.00, 0.06, 0.33, 0.45, -0.025,  -0.020),
    ("B2", 0.00, 0.06, 0.33, 0.45,  0.020,   0.025),
]


def solve(args):
    label, N, T, R, h1, h2, p1, p2, q1, q2, threads = args
    L = 2.0 / N
    t0 = time.time()
    try:
        prob = WhiteConvexProblem(
            N=N, L=L, T=T, R=R,
            h_1=h1, h_2=h2, p_1=p1, p_2=p2,
            q_1=q1, q_2=q2,
        )
        r = prob.solve_and_verify(
            solver=cp.MOSEK, verbose=False,
            mosek_params={"MSK_IPAR_NUM_THREADS": int(threads)},
        )
        return {
            "label": label, "N": N, "T": T, "R": R,
            "ok": True,
            "omega": r["primal_optimum"],
            "max_violation": r["max_violation"],
            "elapsed_s": time.time() - t0,
        }
    except Exception as e:
        import traceback
        return {
            "label": label, "N": N, "T": T, "R": R,
            "ok": False, "error": f"{type(e).__name__}: {e}",
            "traceback": traceback.format_exc(),
            "elapsed_s": time.time() - t0,
        }


def _line(r):
    if r.get("ok"):
        verdict = "PASS" if r["omega"] >= THRESHOLD else "FAIL"
        print(f"  {r['label']:>4s}  N={r['N']:>5d} T={r['T']:>5d}  "
              f"omega={r['omega']:.6f}  mv={r['max_violation']:.1e}  "
              f"{verdict}  ({r['elapsed_s']:.1f}s)", flush=True)
    else:
        print(f"  {r['label']:>4s}  N={r['N']:>5d} T={r['T']:>5d}  "
              f"ERROR {r.get('error')}", flush=True)


def main():
    out_path = HERE / "table2_bump_results.json"
    all_results = []

    # ---- Batch 1: A6/A7/A9 at N=20000, T=5000 ----
    batch1_labels = {"A6", "A7", "A9"}
    batch1_args = [
        (lbl, 20000, 5000, 10, h1, h2, p1, p2, q1, q2, 8)
        for (lbl, h1, h2, p1, p2, q1, q2) in FAILS
        if lbl in batch1_labels
    ]
    print(f"=== Batch 1: {len(batch1_args)} jobs at N=20000, T=5000  "
          f"(4 workers x 8 threads) ===", flush=True)
    t_b1 = time.time()
    with ProcessPoolExecutor(max_workers=4) as ex:
        futs = {ex.submit(solve, a): a[0] for a in batch1_args}
        for fut in as_completed(futs):
            r = fut.result()
            all_results.append(r)
            _line(r)
    print(f"Batch 1 wall: {time.time()-t_b1:.1f}s", flush=True)

    # ---- Batch 2: ALL 5 at N=50000, T=10000 ----
    batch2_args = [
        (lbl, 50000, 10000, 10, h1, h2, p1, p2, q1, q2, 8)
        for (lbl, h1, h2, p1, p2, q1, q2) in FAILS
    ]
    print(f"\n=== Batch 2: {len(batch2_args)} jobs at N=50000, T=10000  "
          f"(2 workers x 8 threads) ===", flush=True)
    t_b2 = time.time()
    with ProcessPoolExecutor(max_workers=2) as ex:
        futs = {ex.submit(solve, a): a[0] for a in batch2_args}
        for fut in as_completed(futs):
            r = fut.result()
            all_results.append(r)
            _line(r)
    print(f"Batch 2 wall: {time.time()-t_b2:.1f}s", flush=True)

    with open(out_path, "w") as f:
        json.dump({"threshold": THRESHOLD, "results": all_results}, f, indent=1)

    # Summary: best omega per label across bumps + paper N=10k baseline.
    BASELINE = {"A6": 0.377160, "A7": 0.374513, "A9": 0.375936,
                "B1": 0.376660, "B2": 0.377907}
    print("\n=== SUMMARY (per row, best Omega across bump levels) ===")
    by_label = {}
    for r in all_results:
        if not r.get("ok"):
            continue
        by_label.setdefault(r["label"], []).append(r)
    for lbl in ["A6", "A7", "A9", "B1", "B2"]:
        runs = by_label.get(lbl, [])
        runs_str = "  ".join(
            f"N={r['N']//1000}k:{r['omega']:.5f}" for r in
            sorted(runs, key=lambda x: x["N"])
        )
        best = max((r["omega"] for r in runs), default=float("-inf"))
        all_runs_str = f"baseline N=10k:{BASELINE[lbl]:.5f}  {runs_str}"
        verdict = "PASS" if best >= THRESHOLD else "FAIL"
        print(f"  {lbl}: {all_runs_str}   best={best:.5f}  {verdict}")
    print(f"\nsaved {out_path}")


if __name__ == "__main__":
    main()
