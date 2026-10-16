"""Probe the worst R=100 cell (h=[0.0272,0.0304], p=[0,0.04], Ω=0.378623)
with several rescue strategies in parallel. Reports the Ω gain per strategy
so we can pick a winner and apply it to all 44 failing cells.

All strategies use R=100, since we already know R=50→100 was saturated.
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

# Parent: h in [0.0272, 0.0304], p in [0, 0.04]; Ω_parent ≈ 0.378623 at N=10000, T=2000, R=100
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
    h_mid = 0.5 * (H_LO + H_HI)
    p_mid = 0.5 * (P_LO + P_HI)

    # Strategy table: each entry is (strategy_label, list of (sub_label, h1,h2,p1,p2,N,T,R))
    jobs = []

    # Baseline reproduction at R=100 (sanity check)
    jobs.append(("baseline",       H_LO, H_HI, P_LO, P_HI, 10000, 2000, 100))

    # Split-h (R=100)
    jobs.append(("split_h_L",      H_LO, h_mid, P_LO, P_HI, 10000, 2000, 100))
    jobs.append(("split_h_R",      h_mid, H_HI, P_LO, P_HI, 10000, 2000, 100))

    # Split-p (R=100)
    jobs.append(("split_p_L",      H_LO, H_HI, P_LO, p_mid, 10000, 2000, 100))
    jobs.append(("split_p_R",      H_LO, H_HI, p_mid, P_HI, 10000, 2000, 100))

    # Bump N=20000 (R=100)
    jobs.append(("bump_N20k",      H_LO, H_HI, P_LO, P_HI, 20000, 4000, 100))

    # Bump T only (T=4000 at N=10000, R=100)
    jobs.append(("bump_T4k",       H_LO, H_HI, P_LO, P_HI, 10000, 4000, 100))

    # Split-h + bump N=20000
    jobs.append(("Nbump_split_h_L", H_LO, h_mid, P_LO, P_HI, 20000, 4000, 100))
    jobs.append(("Nbump_split_h_R", h_mid, H_HI, P_LO, P_HI, 20000, 4000, 100))

    print(f"=== rescue probe on worst cell ===")
    print(f"  parent: h=[{H_LO},{H_HI}] p=[{P_LO},{P_HI}]  R=100 Ω≈0.378623")
    print(f"  running {len(jobs)} jobs in parallel (4 workers)\n", flush=True)

    results = {}
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=4) as ex:
        futs = {ex.submit(solve_one, *j): j[0] for j in jobs}
        for fut in as_completed(futs):
            r = fut.result()
            results[r["label"]] = r
            if r["ok"]:
                tag = "OK" if r["omega"] >= THRESHOLD else "FAIL"
                print(f"  [{r['label']:<18}] N={r['N']} T={r['T']} R={r['R']}  "
                      f"h={r['h']} p={r['p']}  Ω={r['omega']:.6f}  {tag}  "
                      f"({r['elapsed_s']:.0f}s)", flush=True)
            else:
                print(f"  [{r['label']:<18}] ERROR {r['error']}", flush=True)

    print(f"\n  wall={time.time()-t0:.0f}s")

    # Print a strategy summary
    print(f"\n=== strategy verdicts (need both children ≥ {THRESHOLD}) ===")
    def both(a, b):
        ra, rb = results.get(a), results.get(b)
        if not (ra and rb and ra["ok"] and rb["ok"]):
            return None
        return (ra["omega"], rb["omega"])
    def one(a):
        ra = results.get(a)
        if not (ra and ra["ok"]):
            return None
        return ra["omega"]

    def report(name, omegas):
        if omegas is None:
            print(f"  {name:<24}  incomplete")
            return
        if isinstance(omegas, tuple):
            mn = min(omegas)
            verdict = "SAVES" if mn >= THRESHOLD else "fails"
            print(f"  {name:<24}  min Ω = {mn:.6f}  ({omegas[0]:.6f}, {omegas[1]:.6f})  {verdict}")
        else:
            verdict = "SAVES" if omegas >= THRESHOLD else "fails"
            print(f"  {name:<24}  Ω = {omegas:.6f}  {verdict}")

    report("baseline-R100",          one("baseline"))
    report("split-h",                 both("split_h_L","split_h_R"))
    report("split-p",                 both("split_p_L","split_p_R"))
    report("bump-N=20k",              one("bump_N20k"))
    report("bump-T=4k",               one("bump_T4k"))
    report("N=20k + split-h",         both("Nbump_split_h_L","Nbump_split_h_R"))

    out = HERE / "rescue_worst_results.json"
    with open(out, "w") as f:
        json.dump(results, f, indent=1, default=str)
    print(f"\n  saved {out}")


if __name__ == "__main__":
    main()
