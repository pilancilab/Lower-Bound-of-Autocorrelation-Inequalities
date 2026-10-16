"""v51 sweep at higher K."""

from convex_code import AutocorrLowerBoundV51
import time

if __name__ == "__main__":
    print("v51 production sweep:")
    for N, K in [
        (100, 48),
        (100, 64),
        (200, 64),
        (300, 64),
        (200, 96),
    ]:
        t = time.time()
        try:
            prob = AutocorrLowerBoundV51(N=N, K=K)
            out = prob.solve(solver="MOSEK", verbose=False)
            print(f"N={N:5d}, K={K:3d}: Omega={out.Omega:.6f}  ({time.time()-t:.1f}s)")
        except Exception as e:
            print(f"N={N}, K={K}: failed: {e}")
