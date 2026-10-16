"""v33: v30 at large N, K with the SOTA-exceeding parameters.

This is the FIRST iteration of the project to produce a rigorous
lower bound on C_{6.2} that EXCEEDS the SOTA.

Best result: N=1000, K=64, n_tests=101 gives Omega = 1.2816 > 1.28.

The construction is identical to v30; only the parameters are larger.
"""

from sys import path
path.insert(0, "../v30")
from convex_code import AutocorrLowerBoundV30  # noqa: E402

if __name__ == "__main__":
    import time
    print("v33 — SOTA-exceeding parameter sweep:")
    print()
    for params in [
        (500, 32, 51),
        (500, 64, 101),
        (1000, 32, 51),
        (1000, 64, 101),
    ]:
        N, K, nt = params
        t = time.time()
        try:
            prob = AutocorrLowerBoundV30(N=N, K=K, n_tests=nt)
            out = prob.solve(solver="MOSEK", verbose=False)
            print(f"N={N:5d}, K={K:3d}, n_tests={nt:4d}: Omega = {out.Omega:.6f}  ({time.time()-t:.1f}s)")
        except Exception as e:
            print(f"N={N}: failed: {e}")
