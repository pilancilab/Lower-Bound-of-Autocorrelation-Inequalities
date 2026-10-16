"""v34: parameter sweep of v30 at very large N, K to push past v33."""

from sys import path
path.insert(0, "../v30")
from convex_code import AutocorrLowerBoundV30  # noqa: E402

if __name__ == "__main__":
    import time
    print("v34 — extended parameter sweep:")
    for params in [
        (1000, 64, 101),  # v33's best
        (1000, 96, 101),  # increase K
        (2000, 64, 101),  # increase N
        (2000, 96, 101),
    ]:
        N, K, nt = params
        t = time.time()
        try:
            prob = AutocorrLowerBoundV30(N=N, K=K, n_tests=nt)
            out = prob.solve(solver="MOSEK", verbose=False)
            print(f"N={N:5d}, K={K:3d}, n_tests={nt:4d}: Omega = {out.Omega:.6f}  ({time.time()-t:.1f}s)")
        except Exception as e:
            print(f"N={N}, K={K}: failed: {e}")
