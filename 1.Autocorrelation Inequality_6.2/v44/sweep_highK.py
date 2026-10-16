"""v44 high-K sweep where v30 grid sampling fails."""

from convex_code import AutocorrLowerBoundV44
import time

if __name__ == "__main__":
    print("v44 high-K sweep:")
    for N, K in [
        (200, 128),
        (500, 128),
        (300, 128),
    ]:
        t = time.time()
        try:
            prob = AutocorrLowerBoundV44(N=N, K=K)
            out = prob.solve(solver="MOSEK", verbose=False)
            print(f"N={N:5d}, K={K:3d}: Omega={out.Omega:.6f}  ({time.time()-t:.1f}s)")
        except Exception as e:
            print(f"N={N}, K={K}: failed: {e}")
