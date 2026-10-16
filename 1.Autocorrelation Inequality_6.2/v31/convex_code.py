"""v31: parameter sweep of v30 to push toward / past SOTA 1.28."""

from sys import path
path.insert(0, "../v30")
from convex_code import AutocorrLowerBoundV30  # noqa: E402

if __name__ == "__main__":
    import time
    print("v31 parameter sweep of v30:")
    print()
    for params in [
        (16, 16, 33), (32, 32, 51), (50, 32, 51),
        (100, 32, 51), (200, 32, 51), (500, 32, 51),
        (1000, 32, 51), (2000, 32, 51),
    ]:
        N, K, nt = params
        t = time.time()
        try:
            prob = AutocorrLowerBoundV30(N=N, K=K, n_tests=nt)
            out = prob.solve(solver="MOSEK", verbose=False)
            print(f"N={N:5d}, K={K:3d}, n_tests={nt:3d}: Omega = {out.Omega:.6f}  ({time.time()-t:.1f}s)")
        except Exception as e:
            print(f"N={N}: failed: {e}")
