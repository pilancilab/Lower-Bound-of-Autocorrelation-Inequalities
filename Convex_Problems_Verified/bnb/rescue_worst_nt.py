"""Probe the worst R=100 cell with N/T bumps at R=50 (no splitting).

Parent cell: h=[0.0272,0.0304], p=[0,0.04]  (Ω≈0.378623 at N=10000,T=2000,R=100)
Need Ω >= 0.3795. Splits and R=100 only moved Ω by ~1e-5, so we drop R to 50
and push N, T instead.
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

Q1, Q2 = -1.000, -0.050
THRESHOLD = 0.3795

H_LO, H_HI = 0.0272, 0.0304
P_LO, P_HI = 0.0000, 0.0400


def solve_one(label, h1, h2, p1, p2, N, T, R):
    L = 2.0 / N
    t0 = time.time()
    try:
        prob = WhiteConvexProblem(N=N, L=L, T=T, R=R,
                                   h_1=h1, h_2=h2, p_1=p1, p_2=p2,
                                   q_1=Q1, q_2=Q2)
        r = prob.solve_and_verify(solver=cp.MOSEK, verbose=False,
                                   mosek_params={"MSK_IPAR_NUM_THREADS": 8})
        return {"label": label, "ok": True,
                "h": (h1, h2), "p": (p1, p2), "N": N, "T": T, "R": R,
                "omega": r["primal_optimum"],
                "max_violation": r["max_violation"],
                "elapsed_s": time.time() - t0}
    except Exception as e:
        return {"label": label, "ok": False,
                "h": (h1, h2), "p": (p1, p2), "N": N, "T": T, "R": R,
                "error": str(e), "elapsed_s": time.time() - t0}


def main():
    # (label, h1, h2, p1, p2, N, T, R)
    jobs = [
        ("nt_30k_R50",  H_LO, H_HI, P_LO, P_HI,  30000,  6000, 50),
        ("nt_50k_R50",  H_LO, H_HI, P_LO, P_HI,  50000, 10000, 50),
        ("nt_100k_R50", H_LO, H_HI, P_LO, P_HI, 100000, 20000, 50),
    ]

    print(f"=== rescue probe (N/T bump, R=50) ===")
    print(f"  parent: h=[{H_LO},{H_HI}] p=[{P_LO},{P_HI}]")
    print(f"  baseline Ω≈0.378623 at (10000,2000,100); need >= {THRESHOLD}")
    print(f"  running {len(jobs)} jobs in parallel (3 workers, 8 threads each)\n",
          flush=True)

    results = {}
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=3) as ex:
        futs = {ex.submit(solve_one, *j): j[0] for j in jobs}
        for fut in as_completed(futs):
            r = fut.result()
            results[r["label"]] = r
            if r["ok"]:
                tag = "OK" if r["omega"] >= THRESHOLD else "FAIL"
                print(f"  [{r['label']:<12}] N={r['N']} T={r['T']} R={r['R']}  "
                      f"Ω={r['omega']:.6f}  {tag}  "
                      f"viol={r['max_violation']:.2e}  ({r['elapsed_s']:.0f}s)",
                      flush=True)
            else:
                print(f"  [{r['label']:<12}] ERROR {r['error']}  "
                      f"({r['elapsed_s']:.0f}s)", flush=True)

    print(f"\n  wall={time.time()-t0:.0f}s")

    out = HERE / "rescue_worst_nt_results.json"
    with open(out, "w") as f:
        json.dump(results, f, indent=1, default=str)
    print(f"  saved {out}")


if __name__ == "__main__":
    main()
