"""Try both 'smaller cells' and 'larger problem' on A9's single failing cell.

A9 failing cell:
    h in [0, 0.008], p in [0.3, 0.4], q in [0.025, 0.05]
    Omega at L0=(10k, 2k, 10) was 0.377567, needs >= 0.3795 (gap +0.00193).
"""
from __future__ import annotations
import os
os.environ.setdefault("OMP_NUM_THREADS", "8")
os.environ.setdefault("MKL_NUM_THREADS", "8")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "8")

import sys, time, json
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import cvxpy as cp
from EMOP import WhiteConvexProblem

THRESHOLD = 0.3795
Q1, Q2 = 0.025, 0.050

# Approach 1: 4 sub-cells at L0=(10k, 2k, 10).
SUBCELLS = [
    ("smaller_a", 10000, 2000, 10, 0.000, 0.004, 0.30, 0.35),
    ("smaller_b", 10000, 2000, 10, 0.000, 0.004, 0.35, 0.40),
    ("smaller_c", 10000, 2000, 10, 0.004, 0.008, 0.30, 0.35),
    ("smaller_d", 10000, 2000, 10, 0.004, 0.008, 0.35, 0.40),
]

# Approach 2: original cell at larger (N, T, R).
BIGGER = [
    ("L1_20k_4k",  20000,  4000, 10, 0.000, 0.008, 0.30, 0.40),
    ("L2_50k_10k", 50000, 10000, 10, 0.000, 0.008, 0.30, 0.40),
]


def solve_box(args):
    label, N, T, R, h1, h2, p1, p2, threads = args
    L = 2.0 / N
    t0 = time.time()
    try:
        prob = WhiteConvexProblem(N=N, L=L, T=T, R=R,
                                   h_1=h1, h_2=h2, p_1=p1, p_2=p2,
                                   q_1=Q1, q_2=Q2)
        r = prob.solve_and_verify(solver=cp.MOSEK, verbose=False,
                                   mosek_params={"MSK_IPAR_NUM_THREADS": int(threads)})
        return {"label": label, "N": N, "T": T, "R": R,
                "h1": h1, "h2": h2, "p1": p1, "p2": p2,
                "ok": True,
                "omega": r["primal_optimum"],
                "max_violation": r["max_violation"],
                "elapsed_s": time.time() - t0}
    except Exception as e:
        import traceback
        return {"label": label, "ok": False, "error": str(e),
                "traceback": traceback.format_exc()}


def main():
    # Smaller cells at 8 threads/worker (4 small jobs, 4 workers)
    small_args = [(*r, 8) for r in SUBCELLS]
    # Bigger at 8 threads/worker (2 big jobs, 2 workers)
    big_args = [(*r, 8) for r in BIGGER]
    all_args = small_args + big_args   # mix; run in parallel
    print(f"=== running {len(all_args)} jobs in parallel ===", flush=True)
    t0 = time.time()
    results = []
    with ProcessPoolExecutor(max_workers=6) as ex:
        futs = {ex.submit(solve_box, a): a[0] for a in all_args}
        for fut in as_completed(futs):
            r = fut.result()
            results.append(r)
            if r["ok"]:
                verdict = "PASS" if r["omega"] >= THRESHOLD else "FAIL"
                print(f"  {r['label']:>12s}  N={r['N']:>5d}  T={r['T']:>5d}  "
                      f"h=({r['h1']},{r['h2']})  p=({r['p1']},{r['p2']})  "
                      f"omega={r['omega']:.6f}  mv={r['max_violation']:.1e}  "
                      f"{verdict}  ({r['elapsed_s']:.1f}s)", flush=True)
            else:
                print(f"  {r['label']:>12s}  ERROR {r.get('error')}", flush=True)
    print(f"\nWall: {time.time() - t0:.1f}s")

    # Summary
    print("\n=== SUMMARY ===")
    smaller = [r for r in results if r['label'].startswith('smaller') and r.get('ok')]
    bigger  = [r for r in results if not r['label'].startswith('smaller') and r.get('ok')]
    if smaller:
        all_pass_small = all(r['omega'] >= THRESHOLD for r in smaller)
        omin = min(r['omega'] for r in smaller)
        omax = max(r['omega'] for r in smaller)
        print(f"  SMALLER cells (4 at L0): omega range [{omin:.5f}, {omax:.5f}]  "
              f"-> {'ALL PASS' if all_pass_small else 'SOME FAIL'}")
    for r in bigger:
        verdict = "PASS" if r['omega'] >= THRESHOLD else "FAIL"
        print(f"  BIGGER {r['label']:>12s}: omega={r['omega']:.6f}  -> {verdict}")
    print(f"\n  Threshold: {THRESHOLD}")
    print(f"  Original cell at L0=(10k,2k,10): omega=0.377567 (gap +0.001933)")

    out = Path(__file__).resolve().parent / "fix_a9_results.json"
    with open(out, "w") as f:
        json.dump({"threshold": THRESHOLD, "results": results}, f, indent=1)
    print(f"\nsaved {out}")


if __name__ == "__main__":
    main()
