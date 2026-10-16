"""v48 is just v44 with a larger K. Re-export for clarity."""

from sys import path
path.insert(0, "../v44")
from convex_code import AutocorrLowerBoundV44 as AutocorrLowerBoundV48  # noqa: E402, F401


if __name__ == "__main__":
    import time

    print("v48 (= v44 at K >= 128):")
    for N, K in [
        (500, 128),
        (1000, 128),
        (500, 160),
    ]:
        t = time.time()
        try:
            prob = AutocorrLowerBoundV48(N=N, K=K)
            out = prob.solve(solver="MOSEK", verbose=False)
            print(f"N={N:5d}, K={K:3d}: Omega={out.Omega:.6f}  ({time.time()-t:.1f}s)")
        except Exception as e:
            print(f"N={N}, K={K}: failed: {e}")
