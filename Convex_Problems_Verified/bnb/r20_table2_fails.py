"""Re-run the 5 failing Table 2 rows at N=20000, T=4000, R=20.

Idea: R=20 doubles the Bochner-Herglotz Toeplitz block (was 11x11 at R=10,
now 21x21). N stays at 20k and T at 4k. Tests whether more spectral
modes in the PSD block (rather than more N) closes the gap to 0.38.
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
N, T, R = 20000, 4000, 20

# (label, h1, h2, p1, p2, q1, q2)
FAILS = [
    ("A6", 0.00, 0.08, 0.00, 1.00, -1.000, -0.050),
    ("A7", 0.00, 0.08, 0.00, 1.00, -0.050, -0.025),
    ("A9", 0.00, 0.08, 0.00, 1.00,  0.025,  0.050),
    ("B1", 0.00, 0.06, 0.33, 0.45, -0.025, -0.020),
    ("B2", 0.00, 0.06, 0.33, 0.45,  0.020,  0.025),
]


def solve(args):
    label, h1, h2, p1, p2, q1, q2, threads = args
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
        return {"label": label, "ok": True,
                "omega": r["primal_optimum"],
                "max_violation": r["max_violation"],
                "elapsed_s": time.time() - t0}
    except Exception as e:
        import traceback
        return {"label": label, "ok": False,
                "error": f"{type(e).__name__}: {e}",
                "traceback": traceback.format_exc(),
                "elapsed_s": time.time() - t0}


def main():
    print(f"=== {len(FAILS)} rows at N={N}, T={T}, R={R}  "
          f"(4 workers x 8 threads) ===", flush=True)
    args = [(lbl, h1, h2, p1, p2, q1, q2, 8) for (lbl, h1, h2, p1, p2, q1, q2) in FAILS]
    results = []
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=4) as ex:
        futs = {ex.submit(solve, a): a[0] for a in args}
        for fut in as_completed(futs):
            r = fut.result()
            results.append(r)
            if r.get("ok"):
                verdict = "PASS" if r["omega"] >= THRESHOLD else "FAIL"
                print(f"  {r['label']:>4s}  omega={r['omega']:.6f}  "
                      f"mv={r['max_violation']:.1e}  "
                      f"{verdict}  ({r['elapsed_s']:.1f}s)", flush=True)
            else:
                print(f"  {r['label']:>4s}  ERROR {r.get('error')}", flush=True)
    print(f"Wall: {time.time()-t0:.1f}s", flush=True)

    BASELINE_10K = {"A6": 0.377160, "A7": 0.374513, "A9": 0.375936,
                    "B1": 0.376660, "B2": 0.377907}
    BUMP_20K = {"A6": 0.377777, "A7": 0.375034, "A9": 0.376854,
                "B1": None, "B2": None}

    print("\n=== SUMMARY (baselines vs R=20) ===")
    print(f"{'Row':>4s}  {'paper(N=10k,R=10)':>18s}  {'N=20k,R=10':>11s}  "
          f"{'N=20k,R=20':>11s}  {'Verdict'}")
    for r in sorted(results, key=lambda x: x["label"]):
        if not r.get("ok"):
            continue
        lbl = r["label"]
        b10 = BASELINE_10K.get(lbl, float("nan"))
        b20 = BUMP_20K.get(lbl)
        b20_str = f"{b20:.5f}" if b20 is not None else "  N/A  "
        verdict = "PASS" if r["omega"] >= THRESHOLD else "FAIL"
        print(f"  {lbl:>4s}  {b10:>18.5f}  {b20_str:>11s}  "
              f"{r['omega']:>11.5f}  {verdict}")

    out = HERE / "table2_r20_results.json"
    with open(out, "w") as f:
        json.dump({"threshold": THRESHOLD, "N": N, "T": T, "R": R,
                   "results": results}, f, indent=1)
    print(f"\nsaved {out}")


if __name__ == "__main__":
    main()
