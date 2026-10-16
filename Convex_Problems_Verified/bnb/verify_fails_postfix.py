"""Re-run the 5 failing Table 2 rows after fixing 4 -> 8 in eqs (5.6),(5.7)
of EMOP.py. Same (N, T, R) as the paper to test whether the bug-fix alone
recovers the 0.38 bound.
"""
from __future__ import annotations
import os
os.environ.setdefault("OMP_NUM_THREADS", "8")
os.environ.setdefault("MKL_NUM_THREADS", "8")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "8")

import json, sys, time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT))

import cvxpy as cp
from EMOP import WhiteConvexProblem

THRESHOLD = 0.38

# Use paper's (N, T, R) per row.
ROWS = [
    ("A6", 10000, 4000, 10, 0.00, 0.08, 0.00, 1.00, -1.000, -0.050),
    ("A7", 10000, 4000, 10, 0.00, 0.08, 0.00, 1.00, -0.050, -0.025),
    ("A9", 10000, 4000, 10, 0.00, 0.08, 0.00, 1.00,  0.025,  0.050),
    ("B1", 20000, 5000, 10, 0.00, 0.06, 0.33, 0.45, -0.025, -0.020),
    ("B2", 20000, 5000, 10, 0.00, 0.06, 0.33, 0.45,  0.020,  0.025),
]

BEFORE_FIX = {"A6": 0.377160, "A7": 0.374513, "A9": 0.375936,
              "B1": 0.376660, "B2": 0.377907}


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
        return {"label": label, "N": N, "T": T, "R": R, "ok": True,
                "omega": r["primal_optimum"],
                "max_violation": r["max_violation"],
                "elapsed_s": time.time() - t0}
    except Exception as e:
        import traceback
        return {"label": label, "N": N, "T": T, "R": R, "ok": False,
                "error": f"{type(e).__name__}: {e}",
                "traceback": traceback.format_exc(),
                "elapsed_s": time.time() - t0}


def main():
    args = [(*r, 8) for r in ROWS]
    results = []
    print(f"=== {len(args)} rows post-fix (4 workers x 8 threads) ===", flush=True)
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=4) as ex:
        futs = {ex.submit(solve, a): a[0] for a in args}
        for fut in as_completed(futs):
            r = fut.result()
            results.append(r)
            if r.get("ok"):
                v = "PASS" if r["omega"] >= THRESHOLD else "FAIL"
                before = BEFORE_FIX.get(r["label"], float("nan"))
                delta = r["omega"] - before
                print(f"  {r['label']:>4s}  N={r['N']:>5d} T={r['T']:>5d}  "
                      f"omega={r['omega']:.6f}  "
                      f"(was {before:.6f}, Δ={delta:+.6f})  "
                      f"mv={r['max_violation']:.1e}  {v}  "
                      f"({r['elapsed_s']:.1f}s)", flush=True)
            else:
                print(f"  {r['label']:>4s}  ERROR {r.get('error')}", flush=True)
    print(f"Wall: {time.time()-t0:.1f}s")
    out = HERE / "fails_postfix_results.json"
    with open(out, "w") as f:
        json.dump({"threshold": THRESHOLD, "results": results}, f, indent=1)
    print(f"saved {out}")


if __name__ == "__main__":
    main()
