"""v45 at K=128 — should beat v30's broken K=128 result."""

from convex_code import AutocorrLowerBoundV45
import time

if __name__ == "__main__":
    print("v45 at K=128 (exact Fejer-Riesz handles high K cleanly):")
    for N, K in [
        (200, 128),
        (300, 128),
        (500, 128),
    ]:
        t = time.time()
        try:
            prob = AutocorrLowerBoundV45(N=N, K=K)
            out = prob.solve(solver="MOSEK", verbose=False)
            print(f"N={N:5d}, K={K:3d}: Omega={out.Omega:.6f}  ({time.time()-t:.1f}s)")
        except Exception as e:
            print(f"N={N}, K={K}: failed: {e}")
