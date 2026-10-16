"""Recursively split A6's 12 failing cells until all sub-cells clear 0.3795."""
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
Q1, Q2 = -1.000, -0.050  # A6's q range
N, T, R = 10000, 2000, 10
WORKERS = 8

# A6's 12 failing cells from the initial 10x10 verification.
FAILING = [
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


def split_2x2(box):
    h_lo, h_hi, p_lo, p_hi = box
    h_mid = 0.5 * (h_lo + h_hi)
    p_mid = 0.5 * (p_lo + p_hi)
    return [
        (h_lo, h_mid, p_lo, p_mid),
        (h_lo, h_mid, p_mid, p_hi),
        (h_mid, h_hi, p_lo, p_mid),
        (h_mid, h_hi, p_mid, p_hi),
    ]


def solve_box(args):
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
                "error": str(e)}


def main():
    all_results = []
    internal_splits = []
    pending = []
    for box in FAILING:
        children = split_2x2(box)
        internal_splits.append(box)
        for c in children:
            pending.append(c)

    round_num = 0
    while pending:
        round_num += 1
        print(f"\n=== Round {round_num}: {len(pending)} sub-cells ===", flush=True)
        t0 = time.time()
        results = []
        with ProcessPoolExecutor(max_workers=WORKERS) as ex:
            futs = {ex.submit(solve_box, b): b for b in pending}
            for fut in as_completed(futs):
                r = fut.result()
                results.append(r)
        results.sort(key=lambda r: (r["h1"], r["p1"]))
        new_pending = []
        n_pass = n_fail = 0
        for r in results:
            if not r["ok"]:
                print(f"  ERROR h=[{r['h1']:.4f},{r['h2']:.4f}] p=[{r['p1']:.4f},{r['p2']:.4f}]  {r.get('error')}", flush=True)
                continue
            verdict = "PASS" if r["omega"] >= THRESHOLD else "fail->split"
            print(f"  h=[{r['h1']:.5f},{r['h2']:.5f}] p=[{r['p1']:.5f},{r['p2']:.5f}]  "
                  f"omega={r['omega']:.6f} mv={r['max_violation']:.1e} {verdict} ({r['elapsed_s']:.1f}s)",
                  flush=True)
            if r["omega"] >= THRESHOLD:
                all_results.append({**r, "status": "verified"})
                n_pass += 1
            else:
                box = (r["h1"], r["h2"], r["p1"], r["p2"])
                internal_splits.append(box)
                for child in split_2x2(box):
                    new_pending.append(child)
                n_fail += 1
        wall = time.time() - t0
        print(f"  Round {round_num}: pass={n_pass}, split={n_fail}  ({wall:.1f}s)", flush=True)
        pending = new_pending

    print(f"\n=== DONE ===")
    print(f"  total verified leaves: {len(all_results)}")
    print(f"  total internal splits: {len(internal_splits)}")
    print(f"  total rounds:          {round_num}")

    out = HERE / "fix_a6_results.json"
    with open(out, "w") as f:
        json.dump({"threshold": THRESHOLD,
                   "verified_leaves": all_results,
                   "internal_splits": [{"h_lo": b[0], "h_hi": b[1],
                                         "p_lo": b[2], "p_hi": b[3]}
                                        for b in internal_splits]},
                  f, indent=1)
    print(f"  saved {out}")


if __name__ == "__main__":
    main()
