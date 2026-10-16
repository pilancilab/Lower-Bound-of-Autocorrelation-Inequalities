"""Verify a zoomed-in rectangle of A6 with 10x10 grid:
    h in [0.008, 0.04]  (cells of width 0.0032)
    p in [0.00, 0.40]   (cells of width 0.04)
    q = (-1, -0.05)     (A6's full q range)
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
import numpy as np
from EMOP import WhiteConvexProblem

THRESHOLD = 0.3795
N, T, R = 10000, 2000, 10
WORKERS = 8

LABEL = "A6_zoom"
H_LO, H_HI = 0.008, 0.040
P_LO, P_HI = 0.000, 0.400
Q1, Q2 = -1.000, -0.050


def solve_cell(args):
    h1, h2, p1, p2 = args
    L = 2.0 / N
    t0 = time.time()
    try:
        prob = WhiteConvexProblem(N=N, L=L, T=T, R=R,
                                   h_1=h1, h_2=h2, p_1=p1, p_2=p2,
                                   q_1=Q1, q_2=Q2)
        r = prob.solve_and_verify(solver=cp.MOSEK, verbose=False,
                                   mosek_params={"MSK_IPAR_NUM_THREADS": 8})
        return {"ok": True, "h1": h1, "h2": h2, "p1": p1, "p2": p2,
                "omega": r["primal_optimum"],
                "max_violation": r["max_violation"],
                "elapsed_s": time.time() - t0}
    except Exception as e:
        return {"ok": False, "h1": h1, "h2": h2, "p1": p1, "p2": p2,
                "error": str(e), "elapsed_s": time.time() - t0}


def main():
    import math
    h_edges = np.linspace(H_LO, H_HI, 11)
    p_edges = np.linspace(P_LO, P_HI, 11)
    cells = []
    for i in range(10):
        for j in range(10):
            cells.append((float(h_edges[i]), float(h_edges[i+1]),
                          float(p_edges[j]), float(p_edges[j+1])))

    print(f"=== {LABEL}: 100 cells | N={N} T={T} R={R} | "
          f"h-cell-width={(H_HI-H_LO)/10:.4f}, p-cell-width={(P_HI-P_LO)/10:.4f} | "
          f"q=({Q1}, {Q2}) | workers={WORKERS} ===", flush=True)

    boxes = []
    t0 = time.time()
    n_done = 0
    with ProcessPoolExecutor(max_workers=WORKERS) as ex:
        futs = {ex.submit(solve_cell, c): c for c in cells}
        for fut in as_completed(futs):
            r = fut.result()
            n_done += 1
            if r["ok"]:
                o = r["omega"]
                if o is None:
                    status = "failed"
                    tag = "NO-OMEGA"
                elif not math.isfinite(o):
                    status = "verified"  # inf == vacuously verified
                    tag = "INF"
                elif o >= THRESHOLD:
                    status = "verified"
                    tag = "OK"
                else:
                    status = "failed"
                    tag = "FAIL"
                boxes.append({
                    "h_lo": r["h1"], "h_hi": r["h2"],
                    "p_lo": r["p1"], "p_hi": r["p2"],
                    "level": 0, "status": status,
                    "omega": o, "max_violation": r["max_violation"],
                    "last_solver_time_s": r["elapsed_s"], "error": None,
                })
                o_str = f"{o:.5f}" if (o is not None and math.isfinite(o)) else str(o)
                if tag != "OK":  # only print non-OK to keep noise down
                    print(f"  [{n_done}/100] h=[{r['h1']:.4f},{r['h2']:.4f}] "
                          f"p=[{r['p1']:.4f},{r['p2']:.4f}]  omega={o_str}  "
                          f"{tag}  ({r['elapsed_s']:.1f}s)", flush=True)
            else:
                boxes.append({"h_lo": r["h1"], "h_hi": r["h2"],
                              "p_lo": r["p1"], "p_hi": r["p2"],
                              "level": 0, "status": "failed",
                              "omega": None, "max_violation": None,
                              "last_solver_time_s": None,
                              "error": r["error"]})
                print(f"  [{n_done}/100] h=[{r['h1']:.4f},{r['h2']:.4f}] "
                      f"p=[{r['p1']:.4f},{r['p2']:.4f}]  SOLVER-FAIL  "
                      f"({r['elapsed_s']:.1f}s)", flush=True)

    wall = time.time() - t0
    import collections
    sc = collections.Counter(b["status"] for b in boxes)
    # Detailed breakdown
    ver_inf = sum(1 for b in boxes if b["status"]=="verified" and b["omega"] is not None and not math.isfinite(b["omega"]))
    ver_fin = sum(1 for b in boxes if b["status"]=="verified" and b["omega"] is not None and math.isfinite(b["omega"]))
    fail_thr = sum(1 for b in boxes if b["status"]=="failed" and b["omega"] is not None)
    fail_solver = sum(1 for b in boxes if b["status"]=="failed" and b["omega"] is None)
    print(f"\n=== {LABEL}: wall={wall:.1f}s ===")
    print(f"  verified pass:   {ver_fin}")
    print(f"  infeasible:      {ver_inf}")
    print(f"  threshold-fail:  {fail_thr}")
    print(f"  solver-err:      {fail_solver}")

    state_path = HERE / "table2_rows_state" / f"{LABEL}.json"
    state_path.parent.mkdir(exist_ok=True)
    with open(state_path, "w") as f:
        json.dump({"boxes": boxes}, f, indent=1)
    print(f"  saved {state_path}")


if __name__ == "__main__":
    main()
