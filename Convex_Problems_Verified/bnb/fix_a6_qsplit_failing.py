"""Apply q-split only to A6's 12 failing cells.

For each failing cell (h_lo, h_hi, p_lo, p_hi, q=(-1, -0.05)),
create two q-children:
  q in (-1.000, -0.525)  and  q in (-0.525, -0.050)
and solve each.
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
N, T, R = 10000, 2000, 10
WORKERS = 8

# A6's 12 failing cells (h_lo, h_hi, p_lo, p_hi)
FAILING_CELLS = [
    (0.0080, 0.0160, 0.00, 0.10),
    (0.0080, 0.0160, 0.10, 0.20),
    (0.0160, 0.0240, 0.00, 0.10),
    (0.0160, 0.0240, 0.10, 0.20),
    (0.0160, 0.0240, 0.20, 0.30),
    (0.0240, 0.0320, 0.00, 0.10),
    (0.0240, 0.0320, 0.10, 0.20),
    (0.0240, 0.0320, 0.20, 0.30),
    (0.0240, 0.0320, 0.30, 0.40),
    (0.0320, 0.0400, 0.00, 0.10),
    (0.0320, 0.0400, 0.10, 0.20),
    (0.0320, 0.0400, 0.20, 0.30),
]

# q-split at midpoint -0.525
Q_CHILDREN = [
    (-1.000, -0.525),
    (-0.525, -0.050),
]


def solve_cell(args):
    h1, h2, p1, p2, q1, q2 = args
    L = 2.0 / N
    t0 = time.time()
    try:
        prob = WhiteConvexProblem(N=N, L=L, T=T, R=R,
                                   h_1=h1, h_2=h2, p_1=p1, p_2=p2,
                                   q_1=q1, q_2=q2)
        r = prob.solve_and_verify(solver=cp.MOSEK, verbose=False,
                                   mosek_params={"MSK_IPAR_NUM_THREADS": 8})
        return {"ok": True, "h1": h1, "h2": h2, "p1": p1, "p2": p2,
                "q1": q1, "q2": q2,
                "omega": r["primal_optimum"],
                "max_violation": r["max_violation"],
                "status_str": r["status"],
                "elapsed_s": time.time() - t0}
    except Exception as e:
        return {"ok": False, "h1": h1, "h2": h2, "p1": p1, "p2": p2,
                "q1": q1, "q2": q2,
                "error": str(e), "elapsed_s": time.time() - t0}


def main():
    import math
    # Generate 24 sub-cells
    args_list = []
    for (h1, h2, p1, p2) in FAILING_CELLS:
        for (q1, q2) in Q_CHILDREN:
            args_list.append((h1, h2, p1, p2, q1, q2))
    print(f"=== {len(args_list)} sub-cells (12 cells x 2 q-children) at "
          f"L0=({N}, {T}, {R}) ===", flush=True)

    t0 = time.time()
    results = []
    with ProcessPoolExecutor(max_workers=WORKERS) as ex:
        futs = {ex.submit(solve_cell, a): a for a in args_list}
        for fut in as_completed(futs):
            r = fut.result()
            results.append(r)

    results.sort(key=lambda r: (r["h1"], r["p1"], r["q1"]))

    n_pass = n_thr_fail = n_inf = n_solver_err = 0
    for r in results:
        if not r.get("ok"):
            n_solver_err += 1
            tag = "SOLVER-FAIL"
            extra = r.get("error", "")[:80]
        else:
            o = r["omega"]
            if o is not None and not math.isfinite(o):
                n_inf += 1; tag = "INFEASIBLE"; extra = f"status={r.get('status_str')}"
            elif o is None:
                n_solver_err += 1; tag = "NO-OMEGA"
                extra = f"status={r.get('status_str')}"
            elif o >= THRESHOLD:
                n_pass += 1; tag = "PASS"; extra = ""
            else:
                n_thr_fail += 1; tag = "FAIL"; extra = ""
        omega_str = (f"{r['omega']:.6f}" if r.get("omega") is not None and
                     math.isfinite(r["omega"]) else str(r.get("omega")))
        print(f"  h=[{r['h1']:.4f},{r['h2']:.4f}]  p=[{r['p1']:.4f},{r['p2']:.4f}]  "
              f"q=[{r['q1']:.3f},{r['q2']:.3f}]  omega={omega_str}  {tag}  "
              f"({r['elapsed_s']:.1f}s) {extra}", flush=True)

    wall = time.time() - t0
    print(f"\n=== SUMMARY (wall {wall:.1f}s) ===")
    print(f"  pass:        {n_pass}")
    print(f"  infeasible:  {n_inf}")
    print(f"  threshold-fail: {n_thr_fail}")
    print(f"  solver-err:  {n_solver_err}")

    # Per-original-cell summary
    print("\n=== Per original cell ===")
    by_cell = {}
    for r in results:
        key = (r["h1"], r["h2"], r["p1"], r["p2"])
        by_cell.setdefault(key, []).append(r)
    for key, children in sorted(by_cell.items()):
        statuses = []
        for c in children:
            o = c.get("omega")
            if not c.get("ok") or (o is None and not math.isfinite(o or 0)):
                statuses.append("solver-err" if not c.get("ok") else "no-omega")
            elif o is not None and not math.isfinite(o):
                statuses.append("inf")
            elif o is None:
                statuses.append("no-omega")
            elif o >= THRESHOLD:
                statuses.append(f"OK({o:.5f})")
            else:
                statuses.append(f"FAIL({o:.5f})")
        all_clear = all(s.startswith("OK") or s == "inf" for s in statuses)
        verdict = "ALL CLEAR" if all_clear else "NEEDS MORE"
        print(f"  h={key[0]:.4f}-{key[1]:.4f} p={key[2]:.4f}-{key[3]:.4f}: "
              f"q1-half={statuses[0]}, q2-half={statuses[1]} -> {verdict}")

    out = HERE / "fix_a6_qsplit_results.json"
    with open(out, "w") as f:
        json.dump({"threshold": THRESHOLD, "results": results}, f, indent=1)
    print(f"\nsaved {out}")


if __name__ == "__main__":
    main()
