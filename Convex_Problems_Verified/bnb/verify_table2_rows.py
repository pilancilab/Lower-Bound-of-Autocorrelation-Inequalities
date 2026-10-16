"""Single-pass 10x10 verification of the 5 Table 2 rows.

For each row:
  - 10x10 grid in (h, p), q FIXED at the row's q range.
  - Solve each cell once at (N, T, R) = (20000, 5000, 10) with Toeplitz.
  - No splits, no bumps. Cells either verify (Omega >= 0.3795) or fail.

After this pass we'll have a complete picture of where each row stands
(how many cells pass, how many fail, and by how much) and can decide
next steps from there.

Per-row state files in bnb/table2_rows_state/<label>.json.
Per-row logs in bnb/table2_rows_logs/<label>.log.
"""
from __future__ import annotations

import os
os.environ.setdefault("OMP_NUM_THREADS", "8")
os.environ.setdefault("MKL_NUM_THREADS", "8")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "8")

import argparse
import json
import sys
import time
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Optional

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT))

import cvxpy as cp
from EMOP import WhiteConvexProblem


THRESHOLD = 0.3795
MARGIN = 0.0
FEAS_TOL = 1e-6
INITIAL_GRID = 10  # 10x10

STATE_DIR = HERE / "table2_rows_state"
LOG_DIR = HERE / "table2_rows_logs"

# (label, h_lo, h_hi, p_lo, p_hi, q_lo, q_hi, levels)
# Each row's `levels` is a list of (N, T, R, workers, threads) bump rungs.
# Bump path tries each rung in order; split if last rung still fails.
# Starting size matches the original 101-cell BnB: L0=(10k,2k,10).
LEVELS_DEFAULT = [
    # Single-pass 10x10 verification. No splits, no bumps.
    # Same (N, T, R) as the original 101-cell interior-box verification.
    (10000, 2000, 10, 8, 8),
]

ROWS = [
    {
        "label": "A6", "h_lo": 0.00, "h_hi": 0.08,
        "p_lo": 0.00, "p_hi": 1.00, "q_lo": -1.000, "q_hi": -0.050,
        "levels": LEVELS_DEFAULT,
    },
    {
        "label": "A7", "h_lo": 0.00, "h_hi": 0.08,
        "p_lo": 0.00, "p_hi": 1.00, "q_lo": -0.050, "q_hi": -0.025,
        "levels": LEVELS_DEFAULT,
    },
    {
        "label": "A9", "h_lo": 0.00, "h_hi": 0.08,
        "p_lo": 0.00, "p_hi": 1.00, "q_lo":  0.025, "q_hi":  0.050,
        "levels": LEVELS_DEFAULT,
    },
    {
        "label": "B1", "h_lo": 0.00, "h_hi": 0.06,
        "p_lo": 0.33, "p_hi": 0.45, "q_lo": -0.025, "q_hi": -0.020,
        "levels": LEVELS_DEFAULT,
    },
    {
        "label": "B2", "h_lo": 0.00, "h_hi": 0.06,
        "p_lo": 0.33, "p_hi": 0.45, "q_lo":  0.020, "q_hi":  0.025,
        "levels": LEVELS_DEFAULT,
    },
]

MIN_H_WIDTH = 3e-4
MIN_P_WIDTH = 6e-4


def passes(omega, mv):
    return (omega is not None
            and mv is not None
            and omega >= THRESHOLD + MARGIN
            and mv <= FEAS_TOL)


@dataclass
class Box:
    h_lo: float
    h_hi: float
    p_lo: float
    p_hi: float
    level: int = 0
    status: str = "pending"  # pending | verified | failed | split
    omega: Optional[float] = None
    max_violation: Optional[float] = None
    last_solver_time_s: Optional[float] = None
    error: Optional[str] = None

    def key(self):
        return (f"h[{self.h_lo:.6f},{self.h_hi:.6f}]"
                f"_p[{self.p_lo:.6f},{self.p_hi:.6f}]")

    def to_dict(self):
        return asdict(self)

    @classmethod
    def from_dict(cls, d):
        return cls(**d)


def solve_cell(args):
    h_lo, h_hi, p_lo, p_hi, q_lo, q_hi, N, T, R, threads = args
    L = 2.0 / N
    t0 = time.time()
    try:
        prob = WhiteConvexProblem(
            N=N, L=L, T=T, R=R,
            h_1=h_lo, h_2=h_hi,
            p_1=p_lo, p_2=p_hi,
            q_1=q_lo, q_2=q_hi,
        )
        r = prob.solve_and_verify(
            solver=cp.MOSEK, verbose=False,
            mosek_params={"MSK_IPAR_NUM_THREADS": int(threads)},
        )
        return {"ok": True,
                "omega": r["primal_optimum"],
                "max_violation": r["max_violation"],
                "elapsed_s": time.time() - t0}
    except Exception as e:
        return {"ok": False,
                "error": f"{type(e).__name__}: {e}",
                "traceback": traceback.format_exc(),
                "elapsed_s": time.time() - t0}


def split_box(box: Box, max_level: int) -> Optional[tuple[Box, Box]]:
    h_w = box.h_hi - box.h_lo
    p_w = box.p_hi - box.p_lo
    can_h = h_w >= 2 * MIN_H_WIDTH
    can_p = p_w >= 2 * MIN_P_WIDTH
    if not (can_h or can_p):
        return None
    if can_h and (not can_p or h_w >= p_w):
        mid = 0.5 * (box.h_lo + box.h_hi)
        return (Box(box.h_lo, mid, box.p_lo, box.p_hi, level=box.level),
                Box(mid, box.h_hi, box.p_lo, box.p_hi, level=box.level))
    mid = 0.5 * (box.p_lo + box.p_hi)
    return (Box(box.h_lo, box.h_hi, box.p_lo, mid, level=box.level),
            Box(box.h_lo, box.h_hi, mid, box.p_hi, level=box.level))


def initial_boxes(row):
    h_edges = np.linspace(row["h_lo"], row["h_hi"], INITIAL_GRID + 1)
    p_edges = np.linspace(row["p_lo"], row["p_hi"], INITIAL_GRID + 1)
    out = []
    for i in range(INITIAL_GRID):
        for j in range(INITIAL_GRID):
            out.append(Box(
                h_lo=float(h_edges[i]), h_hi=float(h_edges[i + 1]),
                p_lo=float(p_edges[j]), p_hi=float(p_edges[j + 1]),
                level=0,
            ))
    return out


def save_state(boxes, path):
    tmp = path.with_suffix(path.suffix + ".tmp")
    with open(tmp, "w") as f:
        json.dump({"boxes": [b.to_dict() for b in boxes]}, f, indent=1)
    tmp.replace(path)


def load_state(path):
    with open(path) as f:
        return [Box.from_dict(b) for b in json.load(f)["boxes"]]


def summarize(boxes):
    from collections import Counter
    sc = Counter(b.status for b in boxes)
    pl = Counter(b.level for b in boxes if b.status == "pending")
    return ", ".join(f"{k}={v}" for k, v in sorted(sc.items())) + \
           " | pending-by-level=" + repr(dict(pl))


def run_level(row, boxes, level, state_path, log_f):
    pending = [b for b in boxes if b.status == "pending" and b.level == level]
    if not pending:
        return False
    N, T, R, workers, threads = row["levels"][level]
    msg = (f"\n[{row['label']}] L{level}: {len(pending)} boxes | "
           f"N={N} T={T} R={R} | workers={workers} threads={threads}")
    print(msg, flush=True)
    log_f.write(msg + "\n"); log_f.flush()

    args_list = [
        (b.h_lo, b.h_hi, b.p_lo, b.p_hi, row["q_lo"], row["q_hi"],
         N, T, R, threads)
        for b in pending
    ]
    completed = 0
    max_level = len(row["levels"]) - 1
    try:
        with ProcessPoolExecutor(max_workers=workers) as ex:
            futs = {ex.submit(solve_cell, a): i for i, a in enumerate(args_list)}
            for fut in as_completed(futs):
                i = futs[fut]
                box = pending[i]
                try:
                    res = fut.result()
                except Exception as e:
                    res = {"ok": False, "error": f"future-raised: {e}"}
                completed += 1
                _apply_result(box, res, boxes, max_level)
                save_state(boxes, state_path)
                _log_box(box, res, completed, len(pending), log_f)
    except Exception as e:
        log_f.write(f"  pool-error: {e}\n"); log_f.flush()
        print(f"  pool-error: {e}", flush=True)

    end = f"[{row['label']}] L{level} done. {summarize(boxes)}"
    print(end, flush=True)
    log_f.write(end + "\n"); log_f.flush()
    return True


def _apply_result(box: Box, res: dict, boxes: list, max_level: int):
    """Single-pass: solve the cell, then either verify or mark failed.
    No splits, no bumps. (After this pass, the user decides next steps.)
    """
    if not res.get("ok"):
        box.error = res.get("error")
        box.status = "failed"
        return
    box.omega = res["omega"]
    box.max_violation = res["max_violation"]
    box.last_solver_time_s = res["elapsed_s"]
    if passes(box.omega, box.max_violation):
        box.status = "verified"
    else:
        box.status = "failed"


def _log_box(box: Box, res: dict, completed: int, total: int, log_f):
    if not res.get("ok"):
        line = f"  [{completed}/{total}] {box.key()} ERROR {res.get('error')}"
    else:
        if box.status == "verified": action = "OK"
        elif box.status == "split": action = "split"
        elif box.status == "failed": action = "FAILED-min-width"
        else: action = f"bump->L{box.level}"
        line = (f"  [{completed}/{total}] {box.key()} "
                f"omega={res['omega']:.5f} mv={res['max_violation']:.1e} "
                f"-> {action} ({res['elapsed_s']:.1f}s)")
    print(line, flush=True)
    log_f.write(line + "\n"); log_f.flush()


def verify_row(row):
    STATE_DIR.mkdir(parents=True, exist_ok=True)
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    state_path = STATE_DIR / f"{row['label']}.json"
    log_path = LOG_DIR / f"{row['label']}.log"

    if state_path.exists():
        boxes = load_state(state_path)
        action = "resume"
    else:
        boxes = initial_boxes(row)
        save_state(boxes, state_path)
        action = "start"

    max_level = len(row["levels"]) - 1
    with open(log_path, "a") as log_f:
        log_f.write(f"\n=== {action} {row['label']}: "
                    f"q=({row['q_lo']}, {row['q_hi']}), "
                    f"initial={len(boxes)} boxes ===\n")
        log_f.flush()
        while True:
            progressed = False
            for level in range(max_level + 1):
                if any(b.status == "pending" and b.level == level for b in boxes):
                    run_level(row, boxes, level, state_path, log_f)
                    progressed = True
                    break
            if not progressed:
                break

    leaves = [b for b in boxes if b.status != "split"]
    verified = sum(1 for b in leaves if b.status == "verified")
    failed   = sum(1 for b in leaves if b.status == "failed")
    return {"label": row["label"],
            "total_leaves": len(leaves),
            "verified": verified, "failed": failed,
            "n_internal_split": sum(1 for b in boxes if b.status == "split")}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--row", type=str, default=None,
                        help="Limit to a single row label (A6, A7, A9, B1, B2).")
    args = parser.parse_args()

    rows = ROWS if args.row is None else [r for r in ROWS if r["label"] == args.row]
    if not rows:
        print(f"No such row: {args.row}")
        return 1

    results = []
    for row in rows:
        print(f"\n========== ROW {row['label']} ==========", flush=True)
        r = verify_row(row)
        results.append(r)
        print(f"  -> verified {r['verified']}/{r['total_leaves']}, "
              f"failed {r['failed']}, splits {r['n_internal_split']}",
              flush=True)

    print("\n========== SUMMARY ==========")
    print(f"{'row':>4s}  {'verified':>8s}/{'total':>5s}  "
          f"{'failed':>6s}  {'splits':>6s}")
    for r in results:
        print(f"{r['label']:>4s}  {r['verified']:>8d}/{r['total_leaves']:>5d}  "
              f"{r['failed']:>6d}  {r['n_internal_split']:>6d}")
    return 0 if all(r["failed"] == 0 for r in results) else 1


if __name__ == "__main__":
    sys.exit(main())
