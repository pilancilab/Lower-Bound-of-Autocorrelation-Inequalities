"""Verify the v51 record at (N=200, K=96) with analytical cell bounds."""

from convex_code import AutocorrLowerBoundV51v
import time

if __name__ == "__main__":
    print("v51_verified at the record parameters:")
    for N, K in [
        (200, 96),
    ]:
        t = time.time()
        try:
            prob = AutocorrLowerBoundV51v(N=N, K=K)
            out = prob.solve(solver="MOSEK", verbose=False)
            print(f"N={N:5d}, K={K:3d}: Omega={out.Omega:.8f}  ({time.time()-t:.1f}s)")
        except Exception as e:
            print(f"N={N}, K={K}: failed: {e}")
