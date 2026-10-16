"""v44 at larger N, K=96 (push beyond 1.287)."""

from convex_code import AutocorrLowerBoundV44
import time

if __name__ == "__main__":
    print("v44 large-N sweep:")
    for N, K in [
        (1500, 96),
        (2000, 96),
        (3000, 96),
    ]:
        t = time.time()
        try:
            prob = AutocorrLowerBoundV44(N=N, K=K)
            out = prob.solve(solver="MOSEK", verbose=False)
            print(f"N={N:5d}, K={K:3d}: Omega={out.Omega:.6f}  ({time.time()-t:.1f}s)")
        except Exception as e:
            print(f"N={N}, K={K}: failed: {e}")
