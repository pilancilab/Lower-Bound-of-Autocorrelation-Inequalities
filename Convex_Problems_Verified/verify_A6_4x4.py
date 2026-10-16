"""Verify A6 region (h∈[0,0.08]×p∈[0,1]×q∈[-1,-0.05]) using 4x4
starter at paper Table-2 config (N=10k, T=4k, R=10). Adapt-split for failures."""
from __future__ import annotations
import os, sys, time, json, traceback
import numpy as np
from concurrent.futures import ProcessPoolExecutor, as_completed
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

REGION = ("A6", 10000, 4000, 10, 0.00, 0.08, 0.00, 1.00, -1.000, -0.050)
GRID = 4


def build_cells():
    row, N, T, R, h_lo, h_hi, p_lo, p_hi, q_lo, q_hi = REGION
    h_edges = np.linspace(h_lo, h_hi, GRID + 1)
    p_edges = np.linspace(p_lo, p_hi, GRID + 1)
    cells = []
    for i in range(GRID):
        for j in range(GRID):
            cells.append({
                "row": row, "i": i, "j": j,
                "label": f"{row}_{i}_{j}",
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
    y_path = os.path.join(HERE, "y_witnesses", f"{cell['label']}.npz")
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
            "method": info.get("method", "strict_cone_margin"),
            "elapsed": time.time() - t0,
        }
    except Exception as e:
        return {
            "cell": cell, "ok": False,
            "error": f"{type(e).__name__}: {e}",
            "trace": traceback.format_exc(),
            "elapsed": time.time() - t0,
        }


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--out", default=os.path.join(HERE, "regions_A6_4x4_results.json"))
    args = ap.parse_args()

    cells = build_cells()
    workers = min(args.workers, len(cells))
    row, N, T, R, *_ = REGION
    print(f"=== A6 4x4 — {len(cells)} cells at (N={N}, T={T}, R={R}), workers={workers} ===", flush=True)

    THRESHOLD = 0.37912
    t0 = time.time()
    by_label = {}
    n_cert = n_below = n_inf = n_err = 0
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futs = {ex.submit(verify_one, c): c["label"] for c in cells}
        for fut in as_completed(futs):
            r = fut.result()
            c = r["cell"]
            label = c["label"]
            by_label[label] = r
            if r.get("ok") and r.get("valid"):
                L = r["L"]
                if L == float("inf"):
                    msg = f"[INF ] {label:>6s}  primal_infeasible"
                    n_inf += 1
                elif L >= THRESHOLD:
                    msg = (f"[CERT] {label:>6s}  h=[{c['h1']:.4f},{c['h2']:.4f}] "
                           f"p=[{c['p1']:.4f},{c['p2']:.4f}]  L={L:.7f}")
                    n_cert += 1
                else:
                    msg = (f"[BELOW] {label:>5s}  h=[{c['h1']:.4f},{c['h2']:.4f}] "
                           f"p=[{c['p1']:.4f},{c['p2']:.4f}]  L={L:.7f}  (< {THRESHOLD})")
                    n_below += 1
            else:
                n_err += 1
                err = r.get("error", "cone audit failed")
                msg = f"[ERR ] {label:>6s}  {str(err)[:80]}"
            elapsed = time.time() - t0
            done = n_cert + n_below + n_inf + n_err
            avg = elapsed / max(1, done)
            est_rem = avg * (len(cells) - done)
            print(f"  [{done:2d}/{len(cells)}] {msg}  ({r['elapsed']:.0f}s)  "
                  f"[tot {elapsed:.0f}s, eta {est_rem/60:.0f}m]", flush=True)

    total = time.time() - t0

    out = []
    for c in cells:
        r = by_label.get(c["label"], {})
        out.append({
            "label": c["label"], "row": c["row"], "i": c["i"], "j": c["j"],
            "h1": c["h1"], "h2": c["h2"],
            "p1": c["p1"], "p2": c["p2"],
            "q1": c["q1"], "q2": c["q2"],
            "N": c["N"], "T": c["T"], "R": c["R"],
            "certified_lower_bound": (r.get("L") if r.get("ok") and r.get("valid") else None),
            "valid": bool(r.get("valid", False)),
            "method": r.get("method"),
            "error": r.get("error"),
        })
    with open(args.out, "w") as f:
        json.dump(out, f, indent=2, default=float)

    print(f"\nwrote {args.out}", flush=True)
    print(f"=== SUMMARY: CERT={n_cert}, INF={n_inf}, BELOW={n_below}, ERR={n_err}, "
          f"wall={total:.0f}s ({total/60:.1f}m) ===", flush=True)


if __name__ == "__main__":
    main()
