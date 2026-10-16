"""R-sweep on the worst cell at fixed (N=30000, T=6000).

The worst cell (h=[0.0272,0.0304], p=[0,0.04]) cleared the 0.3795 bar at
R=50 with Ω=0.379646 (margin +0.000146). Find the smallest R that still
suffices — that's the R we should run the rest of the cells at.
"""
from __future__ import annotations
import os
os.environ.setdefault("OMP_NUM_THREADS", "8")
os.environ.setdefault("MKL_NUM_THREADS", "8")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "8")

import sys, time, json
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import cvxpy as cp
from EMOP import WhiteConvexProblem

THRESHOLD = 0.3795
N, T = 30000, 6000
Q1, Q2 = -1.000, -0.050
H_LO, H_HI = 0.0272, 0.0304
P_LO, P_HI = 0.0000, 0.0400


def solve_one(R):
    L = 2.0 / N
    t0 = time.time()
    try:
        prob = WhiteConvexProblem(N=N, L=L, T=T, R=R,
                                   h_1=H_LO, h_2=H_HI, p_1=P_LO, p_2=P_HI,
                                   q_1=Q1, q_2=Q2)
        r = prob.solve_and_verify(solver=cp.MOSEK, verbose=False,
                                   mosek_params={"MSK_IPAR_NUM_THREADS": 8})
        return {"R": R, "ok": True,
                "omega": r["primal_optimum"],
                "max_violation": r["max_violation"],
                "elapsed_s": time.time() - t0}
    except Exception as e:
        return {"R": R, "ok": False, "error": str(e),
                "elapsed_s": time.time() - t0}


def main():
    R_values = [10, 20, 30, 40, 50]

    print(f"=== R-sweep on worst cell at N={N}, T={T} ===")
    print(f"  parent: h=[{H_LO},{H_HI}] p=[{P_LO},{P_HI}]")
    print(f"  need Ω >= {THRESHOLD}")
    print(f"  R values: {R_values}")
    print(f"  running {len(R_values)} jobs in parallel "
          f"({len(R_values)} workers × 8 threads each)\n", flush=True)

    results = {}
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=len(R_values)) as ex:
        futs = {ex.submit(solve_one, R): R for R in R_values}
        for fut in as_completed(futs):
            r = fut.result()
            results[r["R"]] = r
            if r["ok"]:
                tag = "OK" if r["omega"] >= THRESHOLD else "FAIL"
                margin = r["omega"] - THRESHOLD
                print(f"  [R={r['R']:>2d}]  Ω={r['omega']:.6f}  {tag}  "
                      f"margin={margin:+.6f}  viol={r['max_violation']:.2e}  "
                      f"({r['elapsed_s']:.0f}s)", flush=True)
            else:
                print(f"  [R={r['R']:>2d}]  ERROR {r['error']}  "
                      f"({r['elapsed_s']:.0f}s)", flush=True)

    print(f"\n  wall={time.time()-t0:.0f}s")

    print("\n=== summary ===")
    print("  R   Ω           verdict   time(s)")
    for R in R_values:
        r = results.get(R)
        if r and r["ok"]:
            v = "PASS" if r["omega"] >= THRESHOLD else "fail"
            print(f"  {R:>2d}  {r['omega']:.6f}  {v:>7s}   {r['elapsed_s']:.0f}")
        else:
            print(f"  {R:>2d}  ERROR")

    out = HERE / "rescue_worst_R_sweep_results.json"
    with open(out, "w") as f:
        json.dump(results, f, indent=1, default=str)
    print(f"\n  saved {out}")


if __name__ == "__main__":
    main()
