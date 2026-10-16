"""Smallest v51 at K=128 — since convergence in N is fast."""

from convex_code import AutocorrLowerBoundV51
import time

if __name__ == "__main__":
    print("v51 K=128 tiny N:")
    for N, K in [
        (32, 128),
        (50, 128),
    ]:
        t = time.time()
        try:
            prob = AutocorrLowerBoundV51(N=N, K=K)
            out = prob.solve(solver="MOSEK", verbose=False)
            print(f"N={N:5d}, K={K:3d}: Omega={out.Omega:.6f}  ({time.time()-t:.1f}s)")
        except Exception as e:
            print(f"N={N}, K={K}: failed: {e}")
