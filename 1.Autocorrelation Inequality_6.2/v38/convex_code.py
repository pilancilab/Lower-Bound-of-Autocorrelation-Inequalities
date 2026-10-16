"""v38: project final state — sweep with smaller, faster parameters.

Designed for relatively quick verification of the project's
core result: C_{6.2} >= 1.28+ via the v30 (autocorrelation Fejer
+ multi-point Fejer mean) construction.
"""

from sys import path
path.insert(0, "../v30")
from convex_code import AutocorrLowerBoundV30  # noqa: E402

if __name__ == "__main__":
    import time
    print("v38 — final verification sweep:")
    for params in [
        (50, 32, 51),    # quick warmup
        (100, 32, 51),
        (200, 32, 51),
        (500, 32, 51),
        (1000, 32, 51),
    ]:
        N, K, nt = params
        t = time.time()
        try:
            prob = AutocorrLowerBoundV30(N=N, K=K, n_tests=nt)
            out = prob.solve(solver="MOSEK", verbose=False)
            print(f"N={N:5d}, K={K:3d}: Omega = {out.Omega:.6f}  ({time.time()-t:.1f}s)")
        except Exception as e:
            print(f"N={N}: failed: {e}")
