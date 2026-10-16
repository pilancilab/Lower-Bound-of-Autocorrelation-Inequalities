"""v51 at K=96 with larger N to check convergence."""

from convex_code import AutocorrLowerBoundV51
import time

if __name__ == "__main__":
    print("v51 K=96 convergence check:")
    for N, K in [
        (300, 96),
        (500, 96),
    ]:
        t = time.time()
        try:
            prob = AutocorrLowerBoundV51(N=N, K=K)
            out = prob.solve(solver="MOSEK", verbose=False)
            print(f"N={N:5d}, K={K:3d}: Omega={out.Omega:.6f}  ({time.time()-t:.1f}s)")
        except Exception as e:
            print(f"N={N}, K={K}: failed: {e}")
