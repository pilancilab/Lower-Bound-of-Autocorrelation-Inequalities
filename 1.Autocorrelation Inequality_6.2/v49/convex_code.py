"""v49 = v45 at larger K (parameter extension)."""

from sys import path
path.insert(0, "../v45")
from convex_code import AutocorrLowerBoundV45 as AutocorrLowerBoundV49  # noqa: E402, F401


if __name__ == "__main__":
    import time

    print("v49 (= v45 at K>=128):")
    for N, K in [
        (500, 128),
        (1000, 128),
    ]:
        t = time.time()
        try:
            prob = AutocorrLowerBoundV49(N=N, K=K)
            out = prob.solve(solver="MOSEK", verbose=False)
            print(f"N={N:5d}, K={K:3d}: Omega={out.Omega:.6f}  ({time.time()-t:.1f}s)")
        except Exception as e:
            print(f"N={N}, K={K}: failed: {e}")
