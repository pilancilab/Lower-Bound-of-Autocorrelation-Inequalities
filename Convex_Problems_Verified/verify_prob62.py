"""Rigorous dual certification of Autocorrelation 6.2 at (N=200, K=96)."""
from __future__ import annotations
import os, sys, time, json, traceback, importlib.util
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)


def _load_prob62():
    name = "prob6_2"
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, "prob6.2.py"))
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod.AutocorrLowerBoundV51v


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--N", type=int, default=200)
    ap.add_argument("--K", type=int, default=96)
    ap.add_argument("--out", default=os.path.join(HERE, "prob62_N200_K96_results.json"))
    ap.add_argument("--margin", type=float, default=1e-8)
    ap.add_argument("--solver-eps", type=float, default=1e-10)
    args = ap.parse_args()

    AutocorrLowerBoundV51v = _load_prob62()
    from dual_certification_api import certified_lower_bound_autocorr

    def build_autocorr(N, K):
        return AutocorrLowerBoundV51v(N=N, K=K).problem

    y_path = os.path.join(HERE, "y_witnesses", f"prob62_N{args.N}_K{args.K}.npz")
    print(f"=== prob6.2 N={args.N} K={args.K} cert  "
          f"(margin={args.margin:.0e}, solver_eps={args.solver_eps:.0e}) ===", flush=True)
    t0 = time.time()
    try:
        L, valid, info = certified_lower_bound_autocorr(
            N=args.N, K=args.K,
            build_autocorr=build_autocorr,
            eps_A=1e-12, eps_b=1e-12, eps_c=1e-12,
            margin=args.margin,
            solver_eps=args.solver_eps,
            y_save_path=y_path,
        )
        elapsed = time.time() - t0
        rec = {
            "N": args.N, "K": args.K,
            "certified_lower_bound": float(L),
            "valid": bool(valid),
            "method": info.get("method", "strict_cone_margin"),
            "margin": args.margin,
            "solver_eps": args.solver_eps,
            "elapsed_s": elapsed,
        }
        print(f"[CERT] N={args.N} K={args.K}  L={L:.10f}  "
              f"valid={valid}  ({elapsed:.0f}s = {elapsed/60:.1f}m)", flush=True)
    except Exception as e:
        elapsed = time.time() - t0
        rec = {
            "N": args.N, "K": args.K,
            "certified_lower_bound": None,
            "valid": False,
            "error": f"{type(e).__name__}: {e}",
            "trace": traceback.format_exc(),
            "elapsed_s": elapsed,
        }
        print(f"[ERR ] N={args.N} K={args.K}  {type(e).__name__}: {e}  "
              f"({elapsed:.0f}s)", flush=True)

    with open(args.out, "w") as f:
        json.dump(rec, f, indent=2, default=float)
    print(f"\nwrote {args.out}", flush=True)


if __name__ == "__main__":
    main()
