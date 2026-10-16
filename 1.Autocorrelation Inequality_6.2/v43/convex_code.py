"""v43: continue v30 sweep at intermediate K."""

from sys import path
path.insert(0, "../v30")
from convex_code import AutocorrLowerBoundV30  # noqa: E402

if __name__ == "__main__":
    import time
    print("v43:")
    for params in [
        (2000, 96, 101),
        (1500, 80, 101),
    ]:
        N, K, nt = params
        t = time.time()
        try:
            prob = AutocorrLowerBoundV30(N=N, K=K, n_tests=nt)
            out = prob.solve(solver="MOSEK", verbose=False)
            print(f"N={N:5d}, K={K:3d}, nt={nt:4d}: Omega = {out.Omega:.6f}  ({time.time()-t:.1f}s)")
        except Exception as e:
            print(f"N={N}: failed: {e}")
