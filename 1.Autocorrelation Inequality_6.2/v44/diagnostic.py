"""v44 diagnostic: extract optimal p and check actual (f*f) vs SDP bound."""

from convex_code import AutocorrLowerBoundV44
import numpy as np
import time


def extract_and_verify(N: int, K: int) -> None:
    print(f"\n=== v44 diagnostic at N={N}, K={K} ===")
    t0 = time.time()
    prob = AutocorrLowerBoundV44(N=N, K=K)
    out = prob.solve(solver="MOSEK", verbose=False)
    print(f"SDP Omega = {out.Omega:.6f}  (solve {time.time()-t0:.1f}s)")

    # Extract optimal p
    p_opt = prob.p.value
    L_f = 1.0 / (4 * N)

    # Construct piecewise-constant f on [-1/4, 1/4]
    n_grid = 40000
    xs = np.linspace(-0.5, 0.5, n_grid, endpoint=False)
    dx = 1.0 / n_grid

    f_vals = np.zeros(n_grid)
    for j in range(2 * N):
        x_l = -0.25 + j * L_f
        x_r = x_l + L_f
        mask = (xs >= x_l) & (xs < x_r)
        if mask.sum() > 0:
            f_vals[mask] = p_opt[j] / L_f

    mass = f_vals.sum() * dx
    print(f"Reconstructed mass = {mass:.6f} (should be 1)")

    # Compute (f*f)(t) on a dense grid of t
    ts = np.linspace(-0.5, 0.5, 1001)
    F_vals = np.zeros_like(ts)
    for i, t in enumerate(ts):
        # (f*f)(t) = integral f(x) f(t-x) dx
        f_shifted = np.interp((t - xs) % 1.0 - 0.5, xs, f_vals, left=0, right=0)
        # Use masking to handle support in [-1/4, 1/4]
        t_minus_x = t - xs
        valid = (t_minus_x >= -0.25) & (t_minus_x <= 0.25)
        f_shift = np.zeros_like(xs)
        # Use linear interpolation: f_shift[i] = f(t - xs[i])
        idxs = np.searchsorted(xs, t_minus_x) - 1
        idxs = np.clip(idxs, 0, n_grid - 2)
        alpha = (t_minus_x - xs[idxs]) / dx
        alpha = np.clip(alpha, 0, 1)
        f_shift = (1 - alpha) * f_vals[idxs] + alpha * f_vals[idxs + 1]
        f_shift[~valid] = 0
        F_vals[i] = (f_vals * f_shift).sum() * dx

    max_F = F_vals.max()
    argmax_t = ts[F_vals.argmax()]
    print(f"Actual max (f*f) = {max_F:.6f} at t = {argmax_t:.4f}")
    print(f"Ratio actual/SDP = {max_F/out.Omega:.6f}")
    print(f"Gap: actual - SDP = {max_F - out.Omega:+.6f}")


if __name__ == "__main__":
    extract_and_verify(N=100, K=32)
    extract_and_verify(N=200, K=48)
