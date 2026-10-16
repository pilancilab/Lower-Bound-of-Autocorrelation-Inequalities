"""Verify A6 with q split in half: A6_1 = q in (-1, -0.525), A6_2 = q in (-0.525, -0.05).
10x10 in (h, p) per sub-row. No (h,p) splits, no bumps."""
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

ROWS = [
    {"label": "A6_1", "h": (0.00, 0.08), "p": (0.00, 1.00),
     "q": (-1.000, -0.525)},
    {"label": "A6_2", "h": (0.00, 0.08), "p": (0.00, 1.00),
     "q": (-0.525, -0.050)},
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
                "omega": r["primal_optimum"],
                "max_violation": r["max_violation"],
                "elapsed_s": time.time() - t0}
    except Exception as e:
        return {"ok": False, "h1": h1, "h2": h2, "p1": p1, "p2": p2,
                "error": str(e)}


def verify_row(row):
    import numpy as np
    h_lo, h_hi = row["h"]
    p_lo, p_hi = row["p"]
    q1, q2 = row["q"]
    h_edges = np.linspace(h_lo, h_hi, 11)
    p_edges = np.linspace(p_lo, p_hi, 11)
    cells = []
    for i in range(10):
        for j in range(10):
            cells.append((float(h_edges[i]), float(h_edges[i+1]),
                          float(p_edges[j]), float(p_edges[j+1]),
                          q1, q2))
    print(f"\n=== {row['label']}: 100 cells | N={N} T={T} R={R} | "
          f"q=({q1}, {q2}) | workers={WORKERS} ===", flush=True)
    boxes = []
    t0 = time.time()
    n_done = 0
    with ProcessPoolExecutor(max_workers=WORKERS) as ex:
        futs = {ex.submit(solve_cell, c): c for c in cells}
        for fut in as_completed(futs):
            r = fut.result()
            n_done += 1
            if r["ok"]:
                status = "verified" if r["omega"] >= THRESHOLD else "failed"
                boxes.append({
                    "h_lo": r["h1"], "h_hi": r["h2"],
                    "p_lo": r["p1"], "p_hi": r["p2"],
                    "level": 0, "status": status,
                    "omega": r["omega"], "max_violation": r["max_violation"],
                    "last_solver_time_s": r["elapsed_s"], "error": None,
                })
                tag = "OK" if status == "verified" else "FAIL"
                if status == "failed" or r["omega"] is None or r["omega"] < 0.4:
                    print(f"  [{n_done}/100] h=[{r['h1']:.4f},{r['h2']:.4f}] "
                          f"p=[{r['p1']:.4f},{r['p2']:.4f}] "
                          f"omega={r['omega']:.5f} {tag} ({r['elapsed_s']:.1f}s)",
                          flush=True)
            else:
                boxes.append({"h_lo": r["h1"], "h_hi": r["h2"],
                              "p_lo": r["p1"], "p_hi": r["p2"],
                              "level": 0, "status": "failed",
                              "omega": None, "max_violation": None,
                              "last_solver_time_s": None,
                              "error": r["error"]})
                print(f"  [{n_done}/100] {r['h1']},{r['h2']},{r['p1']},{r['p2']}  "
                      f"ERROR {r['error']}", flush=True)
    wall = time.time() - t0
    import collections
    sc = collections.Counter(b["status"] for b in boxes)
    print(f"\n  {row['label']}: {dict(sc)}  wall={wall:.1f}s")
    state_path = HERE / "table2_rows_state" / f"{row['label']}.json"
    state_path.parent.mkdir(exist_ok=True)
    with open(state_path, "w") as f:
        json.dump({"boxes": boxes}, f, indent=1)
    print(f"  saved {state_path}")
    return boxes


def main():
    for row in ROWS:
        verify_row(row)


if __name__ == "__main__":
    main()
