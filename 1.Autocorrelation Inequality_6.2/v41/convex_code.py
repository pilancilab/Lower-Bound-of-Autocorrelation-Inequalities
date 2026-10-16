"""v41: K=256 sweep per v33 verdict suggestion."""

from sys import path
path.insert(0, "../v30")
from convex_code import AutocorrLowerBoundV30  # noqa: E402

if __name__ == "__main__":
    import time
    print("v41 K=256 sweep:")
    for params in [
        (250, 256, 129),
        (500, 256, 129),
        (1000, 256, 129),
    ]:
        N, K, nt = params
        t = time.time()
        try:
            prob = AutocorrLowerBoundV30(N=N, K=K, n_tests=nt)
            out = prob.solve(solver="MOSEK", verbose=False)
            print(f"N={N:5d}, K={K:3d}, nt={nt:4d}: Omega = {out.Omega:.6f}  ({time.time()-t:.1f}s)")
        except Exception as e:
            print(f"N={N}: failed: {e}")
