"""Compare v30 and v44 at the same small (N, K)."""

from sys import path
path.insert(0, "../v30")
from convex_code import AutocorrLowerBoundV30  # noqa: E402

import time

print("v30 baseline (n_tests=101) at v44 sanity points:")
for N, K in [(8, 8), (16, 16), (16, 32), (24, 24), (32, 32), (50, 32)]:
    t = time.time()
    try:
        prob = AutocorrLowerBoundV30(N=N, K=K, n_tests=101)
        out = prob.solve(solver="MOSEK", verbose=False)
        print(f"N={N:3d}, K={K:3d}: Omega={out.Omega:.6f}  ({time.time()-t:.1f}s)")
    except Exception as e:
        print(f"N={N}, K={K}: failed: {e}")
