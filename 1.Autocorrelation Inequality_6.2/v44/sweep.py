"""v44 production sweep at parameters where v30 is best."""

from convex_code import AutocorrLowerBoundV44
import time

if __name__ == "__main__":
    print("v44 production sweep:")
    for N, K in [
        (100, 64),
        (300, 64),
        (500, 64),
        (1000, 64),
        (500, 96),
        (1000, 96),
    ]:
        t = time.time()
        try:
            prob = AutocorrLowerBoundV44(N=N, K=K)
            out = prob.solve(solver="MOSEK", verbose=False)
            print(f"N={N:5d}, K={K:3d}: Omega={out.Omega:.6f}  ({time.time()-t:.1f}s)")
        except Exception as e:
            print(f"N={N}, K={K}: failed: {e}")
