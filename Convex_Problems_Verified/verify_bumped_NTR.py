"""For BELOW-threshold cells from a prior round, re-verify at bumped (N, T, R).

Reads a results JSON, filters cells whose certified_lower_bound is below
the threshold, re-runs them at the bumped config, writes a new JSON with
N, T, R included.
"""
from __future__ import annotations
import os, sys, time, json, traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

THRESHOLD = 0.37912


def _label_of(rec):
    return rec.get("sub_label") or rec.get("label") or "(no-label)"


def is_below(r):
    L = r.get("certified_lower_bound")
    return L is not None and L < THRESHOLD


def verify_one(cell):
    sys.path.insert(0, HERE)
    from EMOP import WhiteConvexProblem
    from dual_certification_api import certified_lower_bound_emop

    def build_emop(N, L, T, R, h_1, h_2, p_1, p_2, q_1, q_2):
        return WhiteConvexProblem(N, L, T, R, h_1, h_2, p_1, p_2, q_1, q_2).problem

    N, T, R = cell["N"], cell["T"], cell["R"]
    L_step = 2.0 / N
    y_path = os.path.join(HERE, "y_witnesses", f"{cell['bump_label']}.npz")
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
            "residual_ub": float(info.get("residual_float_ub", float("nan"))),
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
    ap.add_argument("--input",  required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--N", type=int, default=20000)
    ap.add_argument("--T", type=int, default=5000)
    ap.add_argument("--R", type=int, default=10)
    args = ap.parse_args()

    prev = json.load(open(args.input))
    below = [r for r in prev if is_below(r)]

    cells = []
    for r in below:
        c = dict(r)
        c["N"] = args.N
        c["T"] = args.T
        c["R"] = args.R
        c["bump_label"] = f"{_label_of(r)}#N{args.N}T{args.T}R{args.R}"
        c["prior_label"] = _label_of(r)
        c["prior_L"] = r["certified_lower_bound"]
        cells.append(c)

    print(f"=== Bump (N, T, R) -> ({args.N}, {args.T}, {args.R}) "
          f"for {len(cells)} BELOW cells ===", flush=True)
    for c in cells:
        print(f"  {c['prior_label']:>20s} (prior L={c['prior_L']:.5f})  "
              f"h=[{c['h1']:.4f},{c['h2']:.4f}] p=[{c['p1']:.4f},{c['p2']:.4f}] "
              f"q=[{c['q1']:+.3f},{c['q2']:+.3f}]", flush=True)

    workers = min(args.workers, len(cells))
    t0 = time.time()
    results = []
    n_cert = n_below = n_inf = n_err = 0
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futs = {ex.submit(verify_one, c): c["bump_label"] for c in cells}
        for fut in as_completed(futs):
            r = fut.result()
            results.append(r)
            c = r["cell"]
            lab = c["bump_label"]
            if r.get("ok") and r.get("valid"):
                L = r["L"]
                if L == float("inf"):
                    msg = f"[INF ] {lab:>30s}  primal_infeasible"
                    n_inf += 1
                elif L >= THRESHOLD:
                    msg = f"[CERT] {lab:>30s}  L={L:.7f}  (was {c['prior_L']:.7f})"
                    n_cert += 1
                else:
                    msg = f"[BELOW] {lab:>29s}  L={L:.7f}  (was {c['prior_L']:.7f})"
                    n_below += 1
            else:
                n_err += 1
                err = r.get("error", "cone audit failed")
                msg = f"[ERR ] {lab:>30s}  {str(err)[:80]}"
            done = n_cert + n_below + n_inf + n_err
            elapsed = time.time() - t0
            avg = elapsed / max(1, done)
            est_rem = avg * (len(cells) - done)
            print(f"  [{done:2d}/{len(cells)}] {msg}  ({r['elapsed']:.0f}s)  "
                  f"[tot {elapsed:.0f}s, eta {est_rem/60:.1f}m]", flush=True)

    total = time.time() - t0

    out = []
    for r in results:
        c = r["cell"]
        out.append({
            "bump_label":   c["bump_label"],
            "prior_label":  c["prior_label"],
            "prior_certified_lower_bound": c["prior_L"],
            "h1": c["h1"], "h2": c["h2"],
            "p1": c["p1"], "p2": c["p2"],
            "q1": c["q1"], "q2": c["q2"],
            "N":  c["N"],  "T":  c["T"],  "R":  c["R"],
            "certified_lower_bound": (r.get("L") if r.get("ok") and r.get("valid")
                                       else None),
            "valid":   bool(r.get("valid", False)),
            "method":  r.get("method"),
            "error":   r.get("error"),
            "elapsed_s": r.get("elapsed"),
        })
    with open(args.output, "w") as f:
        json.dump(out, f, indent=2, default=float)

    print(f"\nwrote {args.output}", flush=True)
    print(f"=== SUMMARY: CERT={n_cert}, INF={n_inf}, BELOW={n_below}, ERR={n_err}, "
          f"wall={total:.0f}s ({total/60:.1f}m) ===", flush=True)


if __name__ == "__main__":
    main()
