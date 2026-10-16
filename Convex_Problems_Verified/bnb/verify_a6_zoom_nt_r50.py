"""Re-solve A6_zoom's 44 R=100-failing cells at (N=30000, T=6000, R=50).

This config saved the worst cell (h=[0.0272,0.0304], p=[0,0.04]) with
Ω=0.379646 — first config to clear the 0.3795 bar. Apply to the rest.
"""
from __future__ import annotations
import os
os.environ.setdefault("OMP_NUM_THREADS", "8")
os.environ.setdefault("MKL_NUM_THREADS", "8")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "8")

import sys, time, json, math, argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import cvxpy as cp
from EMOP import WhiteConvexProblem

THRESHOLD = 0.3795
N, T = 30000, 6000
R = 50  # default; overridden by --R
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
    parser.add_argument("--workers", type=int, default=4,
                        help="Workers in parallel. ~1.5-3 GB RSS each depending on R.")
    parser.add_argument("--R", type=int, default=50,
                        help="R value (number of harmonic constraints). "
                             "Pick smallest that clears the worst cell.")
    parser.add_argument("--smoke", action="store_true",
                        help="Run only first cell, for a sanity check.")
    parser.add_argument("--source", type=str, default="A6_zoom_R100.json",
                        help="State file in table2_rows_state/ to read "
                             "failing cells from.")
    args = parser.parse_args()
    global R
    R = args.R

    src_path = HERE / "table2_rows_state" / args.source
    src = json.load(open(src_path))
    failing = [b for b in src["boxes"]
               if b["status"] == "failed" and b["omega"] is not None]
    already_ok = [b for b in src["boxes"] if b["status"] == "verified"]
    print(f"=== A6_zoom_NT_R50: re-solving {len(failing)} R=100-failing cells "
          f"at (N={N}, T={T}, R={R}) | workers={args.workers} ===", flush=True)
    print(f"  (keeping {len(already_ok)} cells already verified at R=10/50/100)",
          flush=True)

    cells = [(b["h_lo"], b["h_hi"], b["p_lo"], b["p_hi"]) for b in failing]
    if args.smoke:
        cells = cells[:1]
        print(f"  SMOKE: solving just first cell", flush=True)

    new_boxes = list(already_ok)
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
                new_boxes.append({"h_lo": r["h1"], "h_hi": r["h2"],
                                  "p_lo": r["p1"], "p_hi": r["p2"],
                                  "level": 0, "status": status,
                                  "omega": o,
                                  "max_violation": r["max_violation"],
                                  "last_solver_time_s": r["elapsed_s"],
                                  "error": None})
                o_str = f"{o:.5f}" if (o is not None and math.isfinite(o)) else str(o)
                print(f"  [{n_done}/{len(cells)}] h=[{r['h1']:.4f},{r['h2']:.4f}] "
                      f"p=[{r['p1']:.4f},{r['p2']:.4f}]  omega={o_str}  "
                      f"{tag}  ({r['elapsed_s']:.1f}s)", flush=True)
            else:
                new_boxes.append({"h_lo": r["h1"], "h_hi": r["h2"],
                                  "p_lo": r["p1"], "p_hi": r["p2"],
                                  "level": 0, "status": "failed",
                                  "omega": None, "max_violation": None,
                                  "last_solver_time_s": None,
                                  "error": r["error"]})
                print(f"  [{n_done}/{len(cells)}] h=[{r['h1']:.4f},{r['h2']:.4f}] "
                      f"p=[{r['p1']:.4f},{r['p2']:.4f}]  SOLVER-FAIL  "
                      f"({r['elapsed_s']:.1f}s)", flush=True)

    wall = time.time() - t0
    ver_fin = sum(1 for b in new_boxes if b["status"]=="verified" and b["omega"] is not None and math.isfinite(b["omega"]))
    ver_inf = sum(1 for b in new_boxes if b["status"]=="verified" and b["omega"] is not None and not math.isfinite(b["omega"]))
    fail_thr = sum(1 for b in new_boxes if b["status"]=="failed" and b["omega"] is not None)
    fail_solver = sum(1 for b in new_boxes if b["status"]=="failed" and b["omega"] is None)
    print(f"\n=== A6_zoom_NT_R50: wall={wall:.1f}s ===")
    print(f"  verified pass:  {ver_fin}")
    print(f"  infeasible:     {ver_inf}")
    print(f"  threshold-fail: {fail_thr}")
    print(f"  solver-err:     {fail_solver}")

    dest_path = HERE / "table2_rows_state" / f"A6_zoom_NT_R{R}.json"
    with open(dest_path, "w") as f:
        json.dump({"boxes": new_boxes}, f, indent=1)
    print(f"  saved {dest_path}")


if __name__ == "__main__":
    main()
