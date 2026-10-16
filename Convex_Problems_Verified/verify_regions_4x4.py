"""Verify B1, B2, A7, A9 by splitting each into 4x4 sub-cells. Hardcoded.
Uses paper's (N, T, R) per row. Outputs per-sub-cell JSON."""
from __future__ import annotations
import os, sys, time, json, traceback
import numpy as np
from concurrent.futures import ProcessPoolExecutor, as_completed
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)


# (row, N, T, R, h_lo, h_hi, p_lo, p_hi, q_lo, q_hi) — paper Table 2 defs
REGIONS = [
    ("A7", 10000, 4000, 10, 0.00, 0.08, 0.00, 1.00, -0.050, -0.025),
    ("A9", 10000, 4000, 10, 0.00, 0.08, 0.00, 1.00,  0.025,  0.050),
    ("B1", 20000, 5000, 10, 0.00, 0.06, 0.33, 0.45, -0.025, -0.020),
    ("B2", 20000, 5000, 10, 0.00, 0.06, 0.33, 0.45,  0.020,  0.025),
]

GRID = 4  # 4x4 split


def build_cells():
    cells = []
    for row, N, T, R, h_lo, h_hi, p_lo, p_hi, q_lo, q_hi in REGIONS:
        h_edges = np.linspace(h_lo, h_hi, GRID + 1)
        p_edges = np.linspace(p_lo, p_hi, GRID + 1)
        for i in range(GRID):
            for j in range(GRID):
                cells.append({
                    "row": row, "i": i, "j": j,
                    "N": N, "T": T, "R": R,
                    "h1": float(h_edges[i]),     "h2": float(h_edges[i + 1]),
                    "p1": float(p_edges[j]),     "p2": float(p_edges[j + 1]),
                    "q1": q_lo,                  "q2": q_hi,
                })
    return cells


def verify_one(cell):
    sys.path.insert(0, HERE)
    from EMOP import WhiteConvexProblem
    from dual_certification_api import certified_lower_bound_emop

    def build_emop(N, L, T, R, h_1, h_2, p_1, p_2, q_1, q_2):
        return WhiteConvexProblem(N, L, T, R, h_1, h_2, p_1, p_2, q_1, q_2).problem

    N, T, R = cell["N"], cell["T"], cell["R"]
    L_step = 2.0 / N
    label = f"{cell['row']}_{cell['i']}_{cell['j']}"
    y_path = os.path.join(HERE, "y_witnesses", f"{label}.npz")
    t0 = time.time()
    try:
        L, valid, info = certified_lower_bound_emop(
            N=N, L=L_step, T=T, R=R,
            h_1=cell["h1"], h_2=cell["h2"],
            p_1=cell["p1"], p_2=cell["p2"],
            q_1=cell["q1"], q_2=cell["q2"],
            build_emop=build_emop,
            eps_A=1e-12, eps_b=1e-12, eps_c=1e-12,
            margin=1e-8,
            y_save_path=y_path,
        )
        return {
            "cell": cell, "ok": True, "valid": bool(valid),
            "L": float(L),
            "X": float(info.get("X_used", float("nan"))),
            "residual_ub": float(info.get("residual_float_ub", float("nan"))),
            "y_l1_ub": float(info.get("y_l1_ub", float("nan"))),
            "elapsed": time.time() - t0,
        }
    except Exception as e:
        return {
            "cell": cell, "ok": False,
            "error": f"{type(e).__name__}: {e}",
            "trace": traceback.format_exc(),
            "elapsed": time.time() - t0,
        }


def cell_label(c):
    return f"{c['row']}_{c['i']}_{c['j']}"


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--out", default=os.path.join(HERE, "regions_4x4_results.json"))
    args = ap.parse_args()

    cells = build_cells()
    workers = min(args.workers, len(cells))
    print(f"=== {len(cells)} sub-cells across B1, B2, A7, A9  (4x4 each, workers={workers}) ===", flush=True)
    for row, N, T, R, *_ in REGIONS:
        print(f"    {row}: N={N}, T={T}, R={R}", flush=True)

    t0 = time.time()
    by_label = {}
    n_cert = n_fail = 0
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futs = {ex.submit(verify_one, c): cell_label(c) for c in cells}
        for fut in as_completed(futs):
            r = fut.result()
            c = r["cell"]
            label = cell_label(c)
            by_label[label] = r
            if r.get("ok") and r.get("valid"):
                n_cert += 1
                msg = (f"[CERT] {label:>8s}  h=[{c['h1']:.4f},{c['h2']:.4f}] "
                       f"p=[{c['p1']:.4f},{c['p2']:.4f}] q=[{c['q1']:+.3f},{c['q2']:+.3f}]  "
                       f"L={r['L']:.7f}  ({r['elapsed']:.0f}s)")
            else:
                n_fail += 1
                if r.get("ok"):
                    msg = f"[FAIL] {label:>8s}  cone audit failed  ({r['elapsed']:.0f}s)"
                else:
                    msg = f"[ERR ] {label:>8s}  {r.get('error','?')[:90]}  ({r['elapsed']:.0f}s)"
            elapsed = time.time() - t0
            done = n_cert + n_fail
            avg = elapsed / max(1, done)
            est_rem = avg * (len(cells) - done)
            print(f"  [{done:3d}/{len(cells)}] {msg}  [tot {elapsed:.0f}s, eta {est_rem/60:.0f}m]", flush=True)

    total = time.time() - t0

    out = []
    for c in cells:
        r = by_label.get(cell_label(c), {})
        out.append({
            "label": cell_label(c),
            "row": c["row"], "i": c["i"], "j": c["j"],
            "h1": c["h1"], "h2": c["h2"],
            "p1": c["p1"], "p2": c["p2"],
            "q1": c["q1"], "q2": c["q2"],
            "N": c["N"], "T": c["T"], "R": c["R"],
            "certified_lower_bound": (r.get("L") if r.get("ok") and r.get("valid")
                                       else None),
            "valid": bool(r.get("valid", False)),
            "error": r.get("error"),
        })
    with open(args.out, "w") as f:
        json.dump(out, f, indent=2)
    print(f"\nwrote {args.out}", flush=True)
    print(f"=== SUMMARY: CERT={n_cert}, FAIL+ERR={n_fail}, wall={total:.0f}s ({total/60:.1f}m) ===", flush=True)

    # Per-region min L
    print("\nPer-region min certified L:", flush=True)
    for row, *_ in REGIONS:
        certs = [o for o in out if o["row"] == row and o["certified_lower_bound"] is not None]
        if certs:
            mn = min(certs, key=lambda o: o["certified_lower_bound"])
            print(f"  {row}: min L = {mn['certified_lower_bound']:.7f}  at "
                  f"h=[{mn['h1']:.4f},{mn['h2']:.4f}] p=[{mn['p1']:.4f},{mn['p2']:.4f}]  "
                  f"({len(certs)}/16 certified)", flush=True)
        else:
            print(f"  {row}: NO certified cells (0/16)", flush=True)


if __name__ == "__main__":
    main()
