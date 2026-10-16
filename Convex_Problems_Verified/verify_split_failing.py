"""For each cell from regions_4x4_results.json that failed to certify cleanly
(below threshold or solver error — NOT vacuously verified primal-infeasible),
split into 2x2 sub-cells and run the new dual cert API."""
from __future__ import annotations
import os, sys, time, json, traceback
import numpy as np
from concurrent.futures import ProcessPoolExecutor, as_completed
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

THRESHOLD = 0.37912
DEFAULT_INPUT  = os.path.join(HERE, "regions_4x4_results.json")
DEFAULT_OUTPUT = os.path.join(HERE, "regions_split_2x2_results.json")


def _label_of(parent):
    """Resolve a parent label whether it came from a 4x4 JSON ('label'),
    a previously-split JSON ('sub_label'), or a bumped-NTR JSON ('bump_label',
    with a '#N...T...R...' suffix that we strip)."""
    lab = parent.get("sub_label") or parent.get("label") or parent.get("bump_label")
    if lab and "#" in lab:
        lab = lab.split("#", 1)[0]
    return lab or "(no-label)"


def split_2x2(parent):
    """Split (h, p) box in half along each axis -> 4 sub-cells."""
    hmid = 0.5 * (parent["h1"] + parent["h2"])
    pmid = 0.5 * (parent["p1"] + parent["p2"])
    out = []
    base_label = _label_of(parent)
    for ii, (h_lo, h_hi) in enumerate([(parent["h1"], hmid), (hmid, parent["h2"])]):
        for jj, (p_lo, p_hi) in enumerate([(parent["p1"], pmid), (pmid, parent["p2"])]):
            sub = dict(parent)
            sub["h1"], sub["h2"] = h_lo, h_hi
            sub["p1"], sub["p2"] = p_lo, p_hi
            sub["parent_label"] = base_label
            sub["sub_label"] = f"{base_label}.{ii}{jj}"
            out.append(sub)
    return out


def is_failing(r):
    """A cell is 'failing' if it didn't certify above threshold but is NOT
    vacuously verified by primal infeasibility (the 'zero-size array' bug
    proxy)."""
    if r.get("certified_lower_bound") is not None:
        return r["certified_lower_bound"] < THRESHOLD
    # No certified bound -> some error
    err = (r.get("error") or "")
    if "zero-size array" in err:
        return False  # vacuously verified (infeasible cell)
    return True       # genuine solver failure


def verify_one(sub):
    sys.path.insert(0, HERE)
    from EMOP import WhiteConvexProblem
    from dual_certification_api import certified_lower_bound_emop

    def build_emop(N, L, T, R, h_1, h_2, p_1, p_2, q_1, q_2):
        return WhiteConvexProblem(N, L, T, R, h_1, h_2, p_1, p_2, q_1, q_2).problem

    N, T, R = sub["N"], sub["T"], sub["R"]
    L_step = 2.0 / N
    y_path = os.path.join(HERE, "y_witnesses", f"{sub['sub_label']}.npz")
    t0 = time.time()
    try:
        L, valid, info = certified_lower_bound_emop(
            N=N, L=L_step, T=T, R=R,
            h_1=sub["h1"], h_2=sub["h2"],
            p_1=sub["p1"], p_2=sub["p2"],
            q_1=sub["q1"], q_2=sub["q2"],
            build_emop=build_emop,
            eps_A=1e-12, eps_b=1e-12, eps_c=1e-12,
            margin=1e-8,
            y_save_path=y_path,
        )
        return {
            "sub": sub, "ok": True, "valid": bool(valid),
            "L": float(L),
            "X": float(info.get("X_used", float("nan"))),
            "residual_ub": float(info.get("residual_float_ub", float("nan"))),
            "method": info.get("method", "strict_cone_margin"),
            "elapsed": time.time() - t0,
        }
    except Exception as e:
        return {
            "sub": sub, "ok": False,
            "error": f"{type(e).__name__}: {e}",
            "trace": traceback.format_exc(),
            "elapsed": time.time() - t0,
        }


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--input",  default=DEFAULT_INPUT,
                    help="Parent results JSON to read (default regions_4x4_results.json).")
    ap.add_argument("--output", default=DEFAULT_OUTPUT,
                    help="Output JSON path.")
    args = ap.parse_args()

    parents = json.load(open(args.input))
    failing = [p for p in parents if is_failing(p)]
    subs = []
    for p in failing:
        subs.extend(split_2x2(p))

    print(f"=== Split 2x2 for {len(failing)} failing parent cells -> {len(subs)} sub-cells ===", flush=True)
    print("Failing parents:", flush=True)
    for p in failing:
        L = p.get("certified_lower_bound")
        L_str = f"L={L:.5f}" if L is not None else "ERR"
        err_kind = "below" if L is not None else ("solver_err" if p.get("error") else "?")
        print(f"  {_label_of(p):>12s}  h=[{p['h1']:.4f},{p['h2']:.4f}] p=[{p['p1']:.4f},{p['p2']:.4f}]  {L_str}  ({err_kind})", flush=True)

    workers = min(args.workers, len(subs))
    t0 = time.time()
    results = []
    n_cert = n_inf = n_below = n_err = 0
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futs = {ex.submit(verify_one, sub): sub["sub_label"] for sub in subs}
        for fut in as_completed(futs):
            r = fut.result()
            results.append(r)
            sub = r["sub"]
            label = sub["sub_label"]
            if r.get("ok") and r.get("valid"):
                L = r["L"]
                if L == float("inf"):
                    msg = f"[INF ] {label:>14s}  primal_infeasible_cell_excluded"
                    n_inf += 1
                elif L >= THRESHOLD:
                    msg = f"[CERT] {label:>14s}  L={L:.7f}"
                    n_cert += 1
                else:
                    msg = f"[BELOW] {label:>13s}  L={L:.7f}  (< {THRESHOLD})"
                    n_below += 1
            else:
                n_err += 1
                err = r.get("error", "cone audit failed")
                msg = f"[ERR ] {label:>14s}  {str(err)[:80]}"
            elapsed = time.time() - t0
            done = n_cert + n_inf + n_below + n_err
            avg = elapsed / max(1, done)
            est_rem = avg * (len(subs) - done)
            print(f"  [{done:3d}/{len(subs)}] {msg}  ({r['elapsed']:.0f}s)  "
                  f"[tot {elapsed:.0f}s, eta {est_rem/60:.0f}m]", flush=True)

    total = time.time() - t0

    out = []
    for r in results:
        s = r["sub"]
        out.append({
            "sub_label": s["sub_label"],
            "parent_label": s["parent_label"],
            "h1": s["h1"], "h2": s["h2"],
            "p1": s["p1"], "p2": s["p2"],
            "q1": s["q1"], "q2": s["q2"],
            "N": s["N"], "T": s["T"], "R": s["R"],
            "certified_lower_bound": (r.get("L") if r.get("ok") and r.get("valid") else None),
            "valid": bool(r.get("valid", False)),
            "method": r.get("method"),
            "error": r.get("error"),
        })
    with open(args.output, "w") as f:
        json.dump(out, f, indent=2, default=float)

    print(f"\nwrote {args.output}", flush=True)
    print(f"=== SUMMARY: CERT={n_cert}, INF={n_inf}, BELOW={n_below}, ERR={n_err}, "
          f"wall={total:.0f}s ({total/60:.1f}m) ===", flush=True)

    # Per-parent rollup
    print("\nPer-parent rollup:", flush=True)
    by_parent = {}
    for r in out:
        by_parent.setdefault(r["parent_label"], []).append(r)
    for parent_label, subs in sorted(by_parent.items()):
        certs = [s for s in subs if s["certified_lower_bound"] is not None
                  and s["certified_lower_bound"] != float("inf")]
        infs  = [s for s in subs if s["certified_lower_bound"] == float("inf")]
        below = [s for s in subs if s["certified_lower_bound"] is not None
                  and s["certified_lower_bound"] != float("inf")
                  and s["certified_lower_bound"] < THRESHOLD]
        errs  = [s for s in subs if s["error"] is not None]
        min_L = min((s["certified_lower_bound"] for s in certs), default=None)
        print(f"  {parent_label:>8s}: {len(certs)+len(infs)}/{len(subs)} OK "
              f"(cert={len(certs)-len(below)}, inf={len(infs)}, below={len(below)}, err={len(errs)})  "
              f"min L = {min_L if min_L is not None else 'N/A'}", flush=True)


if __name__ == "__main__":
    main()
