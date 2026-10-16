"""Split each R=100 threshold-failing cell in half (along h, the more
sensitive axis at A6_zoom scale) and re-solve the children at R=50.

Carries forward the 56 cells already verified at R=10/R=50/R=100; only
splits and re-solves the 44 R=100-failing cells (88 children).
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
N, T, R = 10000, 2000, 50
Q1, Q2 = -1.000, -0.050

# A6_zoom row extents (for the wider-normalized-dim split rule)
ROW_H_RANGE = 0.04 - 0.008   # 0.032
ROW_P_RANGE = 0.40 - 0.00    # 0.40


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


def split_child_cells(box):
    """Split the box along the wider normalized dim; tie -> h."""
    hw, pw = box["h_hi"] - box["h_lo"], box["p_hi"] - box["p_lo"]
    h_norm, p_norm = hw / ROW_H_RANGE, pw / ROW_P_RANGE
    if h_norm >= p_norm:
        mid = 0.5 * (box["h_lo"] + box["h_hi"])
        return [(box["h_lo"], mid, box["p_lo"], box["p_hi"]),
                (mid, box["h_hi"], box["p_lo"], box["p_hi"])]
    mid = 0.5 * (box["p_lo"] + box["p_hi"])
    return [(box["h_lo"], box["h_hi"], box["p_lo"], mid),
            (box["h_lo"], box["h_hi"], mid, box["p_hi"])]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=4,
                        help="MOSEK threads per worker = 8.")
    parser.add_argument("--smoke", action="store_true")
    args = parser.parse_args()

    src_path = HERE / "table2_rows_state" / "A6_zoom_R100.json"
    src = json.load(open(src_path))
    failing = [b for b in src["boxes"]
               if b["status"] == "failed" and b["omega"] is not None]
    already_ok = [b for b in src["boxes"] if b["status"] == "verified"]
    print(f"=== A6_zoom_split_R50: splitting {len(failing)} R=100-failing cells "
          f"-> {2*len(failing)} children at R={R} | workers={args.workers} ===",
          flush=True)
    print(f"  (keeping {len(already_ok)} cells already verified)", flush=True)

    cells = []
    for b in failing:
        cells.extend(split_child_cells(b))
    if args.smoke:
        cells = cells[:1]
        print(f"  SMOKE: solving just first child", flush=True)

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
                    status, tag = "failed", "NO-OMEGA"
                elif not math.isfinite(o):
                    status, tag = "verified", "INF"
                elif o >= THRESHOLD:
                    status, tag = "verified", "OK"
                else:
                    status, tag = "failed", "FAIL"
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
    print(f"\n=== A6_zoom_split_R50: wall={wall:.1f}s ===")
    print(f"  verified pass:  {ver_fin}")
    print(f"  infeasible:     {ver_inf}")
    print(f"  threshold-fail: {fail_thr}")
    print(f"  solver-err:     {fail_solver}")

    dest_path = HERE / "table2_rows_state" / "A6_zoom_split_R50.json"
    with open(dest_path, "w") as f:
        json.dump({"boxes": new_boxes}, f, indent=1)
    print(f"  saved {dest_path}")


if __name__ == "__main__":
    main()
