"""Branch-and-bound verification driver for EMOP.py.

Strategy (per user, 2026-05-22):
  - Initial 10x10 grid over (h_1,h_2) x (p_1,p_2). q fixed at (-0.02, 0.02).
  - Bump-first: solve each box at L0 (N=10000). If Omega < THRESHOLD,
    re-solve same box at L1 (N=20000), then L2 (N=50000).
  - Split-last: only if L2 still fails, split the box (wider normalized
    dim) and re-run children inheriting the parent's level.
  - Threshold: Omega >= 0.3795, margin 0.

Run:
  python3 bnb/driver.py
  python3 bnb/driver.py --resume-state bnb/state.json
"""

from __future__ import annotations

# Cap BLAS threads in workers (inherited via fork) so 8 MOSEK workers don't
# fight over 64 cores during scipy/numpy compile-time ops.
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
Q_1, Q_2 = -0.02, 0.02

H_LO_INIT, H_HI_INIT = 0.0, 0.06
P_LO_INIT, P_HI_INIT = 0.33, 0.45

INITIAL_GRID_H = 10
INITIAL_GRID_P = 10

# Bump path: (N, T, R) per level.
LEVELS = [
    (10000, 2000, 10),
    (20000, 4000, 10),
    (50000, 10000, 10),
]
MAX_LEVEL = len(LEVELS) - 1

# (workers, mosek_threads_per_worker). MOSEK scales poorly past 8 threads.
LEVEL_CONFIG = [
    (8, 8),
    (4, 8),
    (2, 8),
]

MIN_H_WIDTH = 3e-4
MIN_P_WIDTH = 6e-4

FEAS_TOL = 1e-6

STATE_PATH = HERE / "state.json"
LOG_PATH = HERE / "log.jsonl"


def passes(omega: Optional[float], max_violation: Optional[float]) -> bool:
    if omega is None or max_violation is None:
        return False
    return (omega >= THRESHOLD + MARGIN) and (max_violation <= FEAS_TOL)


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

    def key(self) -> str:
        return (f"h[{self.h_lo:.6f},{self.h_hi:.6f}]"
                f"_p[{self.p_lo:.6f},{self.p_hi:.6f}]")

    def to_dict(self) -> dict:
        return asdict(self)

    @classmethod
    def from_dict(cls, d: dict) -> "Box":
        return cls(**d)


def solve_box(args):
    h_lo, h_hi, p_lo, p_hi, level, threads = args
    N, T, R = LEVELS[level]
    L = 2.0 / N

    t0 = time.time()
    try:
        problem = WhiteConvexProblem(
            N=N, L=L, T=T, R=R,
            h_1=h_lo, h_2=h_hi,
            p_1=p_lo, p_2=p_hi,
            q_1=Q_1, q_2=Q_2,
        )
        mosek_params = {"MSK_IPAR_NUM_THREADS": int(threads)}
        result = problem.solve_and_verify(
            solver=cp.MOSEK,
            verbose=False,
            mosek_params=mosek_params,
        )
        return {
            "ok": True,
            "status": result["status"],
            "omega": result["primal_optimum"],
            "max_violation": result["max_violation"],
            "elapsed_s": time.time() - t0,
        }
    except Exception as e:
        return {
            "ok": False,
            "error": f"{type(e).__name__}: {e}",
            "traceback": traceback.format_exc(),
            "elapsed_s": time.time() - t0,
        }


def split_box(box: Box) -> Optional[tuple[Box, Box]]:
    h_w = box.h_hi - box.h_lo
    p_w = box.p_hi - box.p_lo
    can_split_h = h_w >= 2 * MIN_H_WIDTH
    can_split_p = p_w >= 2 * MIN_P_WIDTH
    if not (can_split_h or can_split_p):
        return None
    h_norm = h_w / 0.06
    p_norm = p_w / 0.12

    if can_split_h and (not can_split_p or h_norm >= p_norm):
        mid = 0.5 * (box.h_lo + box.h_hi)
        return (
            Box(box.h_lo, mid, box.p_lo, box.p_hi, level=box.level),
            Box(mid, box.h_hi, box.p_lo, box.p_hi, level=box.level),
        )
    mid = 0.5 * (box.p_lo + box.p_hi)
    return (
        Box(box.h_lo, box.h_hi, box.p_lo, mid, level=box.level),
        Box(box.h_lo, box.h_hi, mid, box.p_hi, level=box.level),
    )


def initial_boxes() -> list[Box]:
    h_edges = np.linspace(H_LO_INIT, H_HI_INIT, INITIAL_GRID_H + 1)
    p_edges = np.linspace(P_LO_INIT, P_HI_INIT, INITIAL_GRID_P + 1)
    out = []
    for i in range(INITIAL_GRID_H):
        for j in range(INITIAL_GRID_P):
            out.append(Box(
                h_lo=float(h_edges[i]), h_hi=float(h_edges[i + 1]),
                p_lo=float(p_edges[j]), p_hi=float(p_edges[j + 1]),
                level=0,
            ))
    return out


def save_state(boxes: list[Box], path: Path) -> None:
    tmp = path.with_suffix(path.suffix + ".tmp")
    with open(tmp, "w") as f:
        json.dump({"boxes": [b.to_dict() for b in boxes]}, f, indent=1)
    tmp.replace(path)


def load_state(path: Path) -> list[Box]:
    with open(path, "r") as f:
        data = json.load(f)
    return [Box.from_dict(b) for b in data["boxes"]]


def log_event(event: dict) -> None:
    event["ts"] = time.time()
    with open(LOG_PATH, "a") as f:
        f.write(json.dumps(event) + "\n")


def summarize(boxes: list[Box]) -> str:
    from collections import Counter
    status_counts = Counter(b.status for b in boxes)
    pending_by_level = Counter(b.level for b in boxes if b.status == "pending")
    parts = [f"{k}={v}" for k, v in sorted(status_counts.items())]
    return (", ".join(parts)
            + " | pending-by-level=" + dict(pending_by_level).__repr__())


def run_level_iter(boxes: list[Box], level: int, state_path: Path) -> bool:
    pending = [b for b in boxes if b.status == "pending" and b.level == level]
    if not pending:
        return False
    workers, threads = LEVEL_CONFIG[level]
    N, T, R = LEVELS[level]
    print(f"\n=== L{level} iter: {len(pending)} boxes | N={N} T={T} R={R} "
          f"| workers={workers} threads/worker={threads} ===", flush=True)
    log_event({"event": "iter_start", "level": level, "n_boxes": len(pending),
               "N": N, "T": T, "R": R})

    args_list = [
        (b.h_lo, b.h_hi, b.p_lo, b.p_hi, level, threads)
        for b in pending
    ]

    completed = 0
    try:
        with ProcessPoolExecutor(max_workers=workers) as ex:
            futures = {ex.submit(solve_box, a): i
                       for i, a in enumerate(args_list)}
            for fut in as_completed(futures):
                i = futures[fut]
                box = pending[i]
                try:
                    res = fut.result()
                except Exception as e:
                    res = {"ok": False, "error": f"future-raised: {e}"}
                completed += 1
                _apply_result(box, res, boxes, level)
                save_state(boxes, state_path)
                _print_box_line(box, res, completed, len(pending))
    except Exception as e:
        # Pool died (e.g. BrokenProcessPool). Remaining unprocessed boxes
        # keep status=pending and will be retried next outer iteration.
        msg = f"pool-error: {type(e).__name__}: {e}"
        print(f"  !! {msg}", flush=True)
        log_event({"event": "pool_error", "level": level, "error": msg})

    log_event({"event": "iter_done", "level": level, "summary": summarize(boxes)})
    print(f"=== L{level} done. {summarize(boxes)} ===", flush=True)
    return True


def _apply_result(box: Box, res: dict, boxes: list[Box], level: int) -> None:
    if not res.get("ok", False):
        box.error = res.get("error")
        box.status = "failed"
        log_event({"event": "box_error", "box": box.key(), "level": level,
                   "error": box.error,
                   "traceback": res.get("traceback")})
        return
    omega = res["omega"]
    mv = res["max_violation"]
    box.omega = omega
    box.max_violation = mv
    box.last_solver_time_s = res["elapsed_s"]
    if passes(omega, mv):
        box.status = "verified"
        action = "OK"
    else:
        if box.level < MAX_LEVEL:
            box.level += 1
            action = f"bump->L{box.level}"
        else:
            children = split_box(box)
            if children is None:
                box.status = "failed"
                action = "FAILED-min-width"
            else:
                box.status = "split"
                for c in children:
                    boxes.append(c)
                action = f"split->2*L{children[0].level}"
    log_event({"event": "box_done", "box": box.key(), "level": level,
               "omega": omega, "max_violation": mv, "action": action,
               "elapsed_s": res["elapsed_s"]})


def _print_box_line(box: Box, res: dict, completed: int, total: int) -> None:
    if not res.get("ok", False):
        print(f"  [{completed}/{total}] ERROR {box.key()}: {res.get('error')}",
              flush=True)
        return
    omega = res["omega"]
    mv = res["max_violation"]
    elapsed = res["elapsed_s"]
    if box.status == "verified":
        msg = "OK"
    elif box.status == "split":
        msg = "split"
    elif box.status == "failed":
        msg = "FAILED-min-width"
    else:
        msg = f"bump->L{box.level}"
    print(f"  [{completed}/{total}] {box.key()} omega={omega:.5f} "
          f"mv={mv:.1e} -> {msg} ({elapsed:.1f}s)", flush=True)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--resume-state", type=str, default=None)
    parser.add_argument("--max-level", type=int, default=MAX_LEVEL)
    parser.add_argument("--limit-boxes", type=int, default=None,
                        help="For smoke testing: process only first N boxes.")
    args = parser.parse_args()

    state_path = Path(args.resume_state) if args.resume_state else STATE_PATH

    if args.resume_state and state_path.exists():
        boxes = load_state(state_path)
        print(f"Resumed from {state_path}: {len(boxes)} boxes. "
              f"{summarize(boxes)}", flush=True)
    else:
        boxes = initial_boxes()
        if args.limit_boxes is not None:
            boxes = boxes[:args.limit_boxes]
        save_state(boxes, state_path)
        print(f"Fresh start: {len(boxes)} boxes at L0. State={state_path}",
              flush=True)
    log_event({"event": "driver_start", "n_boxes": len(boxes),
               "threshold": THRESHOLD, "max_level": args.max_level})

    while True:
        progressed = False
        for level in range(args.max_level + 1):
            pending_here = any(
                b.status == "pending" and b.level == level for b in boxes
            )
            if pending_here:
                run_level_iter(boxes, level, state_path)
                progressed = True
                break  # restart scan from L0 (splits may add new boxes)
        if not progressed:
            break

    print("\n=== FINAL ===", flush=True)
    print(summarize(boxes), flush=True)
    leaves = [b for b in boxes if b.status != "split"]
    n_verified = sum(1 for b in leaves if b.status == "verified")
    n_failed = sum(1 for b in leaves if b.status == "failed")
    print(f"verified leaves: {n_verified}/{len(leaves)}")
    print(f"failed leaves:   {n_failed}/{len(leaves)}")
    log_event({"event": "driver_done",
               "verified_leaves": n_verified,
               "failed_leaves": n_failed,
               "total_leaves": len(leaves)})
    return 0 if n_failed == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
