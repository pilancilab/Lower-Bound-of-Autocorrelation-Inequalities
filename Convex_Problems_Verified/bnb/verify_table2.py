"""Verify the '0.38' rows of Table 2 (white_paper.txt) using current EMOP.py.

Each row is one convex solve at the paper's (N, T, R). Pass iff Omega >= 0.38.
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

# (label, N, T, R, h1, h2, p1, p2, q1, q2)
ROWS = [
    # N=10000, T=4000, R=10
    ("A1",  10000, 4000, 10, 0.75, 2.00, 0.00, 1.00, -1.000,  1.000),
    ("A2",  10000, 4000, 10, 0.40, 0.75, 0.00, 1.00, -1.000,  1.000),
    ("A3",  10000, 4000, 10, 0.20, 0.40, 0.00, 1.00, -1.000,  1.000),
    ("A4",  10000, 4000, 10, 0.10, 0.20, 0.00, 1.00, -1.000,  1.000),
    ("A5",  10000, 4000, 10, 0.08, 0.10, 0.00, 1.00, -1.000,  1.000),
    ("A6",  10000, 4000, 10, 0.00, 0.08, 0.00, 1.00, -1.000, -0.050),
    ("A7",  10000, 4000, 10, 0.00, 0.08, 0.00, 1.00, -0.050, -0.025),
    ("A8",  10000, 4000, 10, 0.00, 0.08, 0.00, 1.00,  0.050,  1.000),
    ("A9",  10000, 4000, 10, 0.00, 0.08, 0.00, 1.00,  0.025,  0.050),
    ("A10", 10000, 4000, 10, 0.00, 0.08, 0.00, 0.25, -0.025,  0.025),
    ("A11", 10000, 4000, 10, 0.00, 0.08, 0.25, 0.30, -0.025,  0.025),
    ("A12", 10000, 4000, 10, 0.00, 0.08, 0.30, 0.33, -0.025,  0.025),
    ("A13", 10000, 4000, 10, 0.00, 0.08, 0.50, 1.00, -0.025,  0.025),
    ("A14", 10000, 4000, 10, 0.00, 0.08, 0.45, 0.50, -0.025,  0.025),
    ("A15", 10000, 4000, 10, 0.06, 0.08, 0.33, 0.45, -0.025,  0.025),
    # N=20000, T=5000, R=10
    ("B1",  20000, 5000, 10, 0.00, 0.06, 0.33, 0.45, -0.025, -0.020),
    ("B2",  20000, 5000, 10, 0.00, 0.06, 0.33, 0.45,  0.020,  0.025),
]


def solve_row(args):
    label, N, T, R, h1, h2, p1, p2, q1, q2, threads = args
    L = 2.0 / N
    t0 = time.time()
    try:
        prob = WhiteConvexProblem(
            N=N, L=L, T=T, R=R,
            h_1=h1, h_2=h2, p_1=p1, p_2=p2,
            q_1=q1, q_2=q2,
        )
        result = prob.solve_and_verify(
            solver=cp.MOSEK, verbose=False,
            mosek_params={"MSK_IPAR_NUM_THREADS": int(threads)},
        )
        return {
            "label": label,
            "ok": True,
            "N": N, "T": T, "R": R,
            "h1": h1, "h2": h2, "p1": p1, "p2": p2, "q1": q1, "q2": q2,
            "omega": result["primal_optimum"],
            "status": result["status"],
            "max_violation": result["max_violation"],
            "elapsed_s": time.time() - t0,
        }
    except Exception as e:
        import traceback
        return {
            "label": label, "ok": False,
            "N": N, "T": T, "R": R,
            "h1": h1, "h2": h2, "p1": p1, "p2": p2, "q1": q1, "q2": q2,
            "error": f"{type(e).__name__}: {e}",
            "traceback": traceback.format_exc(),
            "elapsed_s": time.time() - t0,
        }


def _print(r):
    if r.get("ok"):
        verdict = "PASS" if r["omega"] >= THRESHOLD else "FAIL"
        print(f"  {r['label']:>4s}  N={r['N']:>5d} T={r['T']:>5d}  "
              f"h=({r['h1']:>5.3f},{r['h2']:>5.3f})  "
              f"p=({r['p1']:>5.3f},{r['p2']:>5.3f})  "
              f"q=({r['q1']:>+6.3f},{r['q2']:>+6.3f})  "
              f"omega={r['omega']:.6f}  mv={r['max_violation']:.1e}  "
              f"{verdict}  ({r['elapsed_s']:.1f}s)", flush=True)
    else:
        print(f"  {r['label']:>4s}  ERROR  {r.get('error')}", flush=True)


def main():
    n10k = [r for r in ROWS if r[1] == 10000]
    n20k = [r for r in ROWS if r[1] == 20000]
    results = []

    print(f"=== {len(n10k)} rows at N=10000, T=4000  "
          f"(workers=8, threads/worker=8) ===", flush=True)
    args = [(*r, 8) for r in n10k]
    with ProcessPoolExecutor(max_workers=8) as ex:
        futs = {ex.submit(solve_row, a): a[0] for a in args}
        for fut in as_completed(futs):
            r = fut.result()
            results.append(r)
            _print(r)

    print(f"\n=== {len(n20k)} rows at N=20000, T=5000  "
          f"(workers=2, threads/worker=8) ===", flush=True)
    args = [(*r, 8) for r in n20k]
    with ProcessPoolExecutor(max_workers=2) as ex:
        futs = {ex.submit(solve_row, a): a[0] for a in args}
        for fut in as_completed(futs):
            r = fut.result()
            results.append(r)
            _print(r)

    out = HERE / "table2_results.json"
    with open(out, "w") as f:
        json.dump({"threshold": THRESHOLD, "results": results}, f, indent=1)

    print("\n=== SUMMARY ===", flush=True)
    passed = [r for r in results if r.get("ok") and r["omega"] >= THRESHOLD]
    failed = [r for r in results if r.get("ok") and r["omega"] <  THRESHOLD]
    errored = [r for r in results if not r.get("ok")]
    print(f"passed:  {len(passed):>2d}/{len(results)}")
    print(f"failed:  {len(failed):>2d}/{len(results)}")
    print(f"errored: {len(errored):>2d}/{len(results)}")
    if failed:
        print("\nfailed rows:")
        for r in sorted(failed, key=lambda x: x["label"]):
            print(f"  {r['label']}: omega={r['omega']:.6f}")
    if errored:
        print("\nerrored rows:")
        for r in errored:
            print(f"  {r['label']}: {r.get('error')}")
    print(f"\nsaved {out}")


if __name__ == "__main__":
    main()
