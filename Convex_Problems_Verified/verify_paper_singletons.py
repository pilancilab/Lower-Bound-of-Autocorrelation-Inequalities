"""Dual verification of the 12 single-cell paper rows that are NOT A6/A7/A9.
Hardcoded cells from white_paper Table 2; calls certified_lower_bound_emop
in parallel via ProcessPoolExecutor. No CSV / JSON reads."""
from __future__ import annotations
import os, sys, time, json, traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)


# (label, N, T, R, h1, h2, p1, p2, q1, q2) — A1-A5, A8, A10-A15
ROWS = [
    ("A1",  10000, 4000, 10, 0.75, 2.00, 0.00, 1.00, -1.000,  1.000),
    ("A2",  10000, 4000, 10, 0.40, 0.75, 0.00, 1.00, -1.000,  1.000),
    ("A3",  10000, 4000, 10, 0.20, 0.40, 0.00, 1.00, -1.000,  1.000),
    ("A4",  10000, 4000, 10, 0.10, 0.20, 0.00, 1.00, -1.000,  1.000),
    ("A5",  10000, 4000, 10, 0.08, 0.10, 0.00, 1.00, -1.000,  1.000),
    ("A8",  10000, 4000, 10, 0.00, 0.08, 0.00, 1.00,  0.050,  1.000),
    ("A10", 10000, 4000, 10, 0.00, 0.08, 0.00, 0.25, -0.025,  0.025),
    ("A11", 10000, 4000, 10, 0.00, 0.08, 0.25, 0.30, -0.025,  0.025),
    ("A12", 10000, 4000, 10, 0.00, 0.08, 0.30, 0.33, -0.025,  0.025),
    ("A13", 10000, 4000, 10, 0.00, 0.08, 0.50, 1.00, -0.025,  0.025),
    ("A14", 10000, 4000, 10, 0.00, 0.08, 0.45, 0.50, -0.025,  0.025),
    ("A15", 10000, 4000, 10, 0.06, 0.08, 0.33, 0.45, -0.025,  0.025),
]


def verify_one(row):
    sys.path.insert(0, HERE)
    from EMOP import WhiteConvexProblem
    from dual_certification_api import certified_lower_bound_emop

    def build_emop(N, L, T, R, h_1, h_2, p_1, p_2, q_1, q_2):
        return WhiteConvexProblem(N, L, T, R, h_1, h_2, p_1, p_2, q_1, q_2).problem

    label, N, T, R, h1, h2, p1, p2, q1, q2 = row
    L_step = 2.0 / N
    y_path = os.path.join(HERE, "y_witnesses", f"singleton_{label}.npz")
    t0 = time.time()
    try:
        L, valid, info = certified_lower_bound_emop(
            N=N, L=L_step, T=T, R=R,
            h_1=h1, h_2=h2, p_1=p1, p_2=p2, q_1=q1, q_2=q2,
            build_emop=build_emop,
            eps_A=1e-12, eps_b=1e-12, eps_c=1e-12,
            margin=1e-8,
            y_save_path=y_path,
        )
        return {
            "label": label,
            "ok": True, "valid": bool(valid),
            "L": float(L), "X": float(info.get("X_used", float("nan"))),
            "residual_ub": float(info.get("residual_float_ub", float("nan"))),
            "y_l1_ub": float(info.get("y_l1_ub", float("nan"))),
            "row": row, "elapsed": time.time() - t0,
        }
    except Exception as e:
        return {
            "label": label, "ok": False,
            "error": f"{type(e).__name__}: {e}",
            "trace": traceback.format_exc(),
            "row": row, "elapsed": time.time() - t0,
        }


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--out", default=os.path.join(HERE, "paper_singletons_results.json"),
                    help="JSON output with per-cell certified bounds.")
    args = ap.parse_args()
    workers = min(args.workers, len(ROWS))
    print(f"=== Dual cert of {len(ROWS)} paper singletons (A1-A5, A8, A10-A15) ===", flush=True)
    print(f"    workers={workers}, all at (N=10000, T=4000, R=10)", flush=True)

    t0 = time.time()
    n_cert = n_fail = 0
    by_label = {}
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futs = {ex.submit(verify_one, r): r[0] for r in ROWS}
        for fut in as_completed(futs):
            r = fut.result()
            label = r["label"]
            by_label[label] = r
            if r.get("ok") and r.get("valid"):
                n_cert += 1
                row = r["row"]
                print(f"  [CERT] {label:>4s} h=[{row[4]:>5.3f},{row[5]:>5.3f}] "
                      f"p=[{row[6]:>5.3f},{row[7]:>5.3f}] q=[{row[8]:>+6.3f},{row[9]:>+6.3f}]  "
                      f"L={r['L']:.7f}  res_ub={r['residual_ub']:.2e}  "
                      f"y_l1_ub={r['y_l1_ub']:.2e}  X={r['X']:.2e}  "
                      f"({r['elapsed']:.1f}s)", flush=True)
            else:
                n_fail += 1
                if r.get("ok"):
                    print(f"  [FAIL] {label:>4s}  (cone audit failed)  ({r['elapsed']:.1f}s)", flush=True)
                else:
                    print(f"  [ERR ] {label:>4s}  {r.get('error','?')[:90]}  ({r['elapsed']:.1f}s)", flush=True)

    total = time.time() - t0

    # Write summary JSON in paper-row order with just the requested fields.
    out = []
    for row in ROWS:
        label, N, T, R, h1, h2, p1, p2, q1, q2 = row
        r = by_label.get(label, {})
        out.append({
            "label": label,
            "h1": h1, "h2": h2,
            "p1": p1, "p2": p2,
            "q1": q1, "q2": q2,
            "certified_lower_bound": (r.get("L") if r.get("ok") and r.get("valid")
                                       else None),
            "valid": bool(r.get("valid", False)),
            "error": r.get("error"),
        })
    with open(args.out, "w") as f:
        json.dump(out, f, indent=2)
    print(f"\nwrote {args.out}", flush=True)
    print(f"=== SUMMARY: CERT={n_cert}, FAIL+ERR={n_fail}, wall={total:.0f}s ===", flush=True)


if __name__ == "__main__":
    main()
