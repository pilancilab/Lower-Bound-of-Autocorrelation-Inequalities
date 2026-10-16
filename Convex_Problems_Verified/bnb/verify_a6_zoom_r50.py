"""Re-run A6_zoom rectangle at R=50 (much larger Toeplitz SDP block)."""
from __future__ import annotations
import os
os.environ.setdefault("OMP_NUM_THREADS", "8")
os.environ.setdefault("MKL_NUM_THREADS", "8")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "8")

import sys, time, json, argparse, math
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import cvxpy as cp
import numpy as np
from EMOP import WhiteConvexProblem

THRESHOLD = 0.3795
N, T, R = 10000, 2000, 50

LABEL = "A6_zoom_R50"
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
    parser = argparse.ArgumentParser()
    parser.add_argument("--smoke", action="store_true",
                        help="Just solve one of the worst cells to estimate per-cell time.")
    parser.add_argument("--workers", type=int, default=4,
                        help="Parallel workers (lower than 8 due to larger SDP).")
    args = parser.parse_args()

    if args.smoke:
        # Pick a representative worst cell from A6_zoom R=10 (Omega was ~0.378)
        smoke_cell = (0.0240, 0.0272, 0.10, 0.14)
        print(f"=== SMOKE TEST: 1 cell at R={R} ===", flush=True)
        print(f"  cell: h=[{smoke_cell[0]:.4f},{smoke_cell[1]:.4f}] "
              f"p=[{smoke_cell[2]:.4f},{smoke_cell[3]:.4f}]  q=({Q1},{Q2})", flush=True)
        r = solve_cell(smoke_cell)
        if r["ok"]:
            o = r["omega"]
            o_str = f"{o:.6f}" if (o is not None and math.isfinite(o)) else str(o)
            print(f"  omega={o_str}  mv={r['max_violation']:.1e}  "
                  f"elapsed={r['elapsed_s']:.1f}s", flush=True)
            # Compare to R=10 baseline
            from collections import Counter
            try:
                d = json.load(open(HERE / "table2_rows_state" / "A6_zoom.json"))['boxes']
                for b in d:
                    if (abs(b['h_lo'] - smoke_cell[0]) < 1e-6 and
                        abs(b['p_lo'] - smoke_cell[2]) < 1e-6):
                        if b['omega'] is not None:
                            print(f"  R=10 baseline for same cell: omega={b['omega']:.6f}")
                            print(f"  delta (R=50 - R=10): "
                                  f"{r['omega'] - b['omega']:+.6f}")
                        break
            except Exception:
                pass
        else:
            print(f"  ERROR: {r['error']}", flush=True)
        return

    h_edges = np.linspace(H_LO, H_HI, 11)
    p_edges = np.linspace(P_LO, P_HI, 11)
    cells = []
    for i in range(10):
        for j in range(10):
            cells.append((float(h_edges[i]), float(h_edges[i+1]),
                          float(p_edges[j]), float(p_edges[j+1])))

    print(f"=== {LABEL}: 100 cells | N={N} T={T} R={R} | "
          f"workers={args.workers} ===", flush=True)

    boxes = []
    t0 = time.time()
    n_done = 0
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        futs = {ex.submit(solve_cell, c): c for c in cells}
        for fut in as_completed(futs):
            r = fut.result()
            n_done += 1
            if r["ok"]:
                o = r["omega"]
                if o is None:
                    status = "failed"; tag = "NO-OMEGA"
                elif not math.isfinite(o):
                    status = "verified"; tag = "INF"
                elif o >= THRESHOLD:
                    status = "verified"; tag = "OK"
                else:
                    status = "failed"; tag = "FAIL"
                boxes.append({"h_lo": r["h1"], "h_hi": r["h2"],
                              "p_lo": r["p1"], "p_hi": r["p2"],
                              "level": 0, "status": status,
                              "omega": o, "max_violation": r["max_violation"],
                              "last_solver_time_s": r["elapsed_s"],
                              "error": None})
                o_str = f"{o:.5f}" if (o is not None and math.isfinite(o)) else str(o)
                if tag != "OK":
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
    from collections import Counter
    sc = Counter(b["status"] for b in boxes)
    ver_inf = sum(1 for b in boxes if b["status"]=="verified" and b["omega"] is not None and not math.isfinite(b["omega"]))
    ver_fin = sum(1 for b in boxes if b["status"]=="verified" and b["omega"] is not None and math.isfinite(b["omega"]))
    fail_thr = sum(1 for b in boxes if b["status"]=="failed" and b["omega"] is not None)
    fail_solver = sum(1 for b in boxes if b["status"]=="failed" and b["omega"] is None)
    print(f"\n=== {LABEL}: wall={wall:.1f}s ===")
    print(f"  verified pass:  {ver_fin}")
    print(f"  infeasible:     {ver_inf}")
    print(f"  threshold-fail: {fail_thr}")
    print(f"  solver-err:     {fail_solver}")

    state_path = HERE / "table2_rows_state" / f"{LABEL}.json"
    state_path.parent.mkdir(exist_ok=True)
    with open(state_path, "w") as f:
        json.dump({"boxes": boxes}, f, indent=1)
    print(f"  saved {state_path}")


if __name__ == "__main__":
    main()
