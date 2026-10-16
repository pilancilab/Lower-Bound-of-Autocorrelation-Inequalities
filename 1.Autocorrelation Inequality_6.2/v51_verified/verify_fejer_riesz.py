"""Sanity-check the Fejer-Riesz constraint derivation.

For a known f (uniform on [-1/4, 1/4]), compute:
  (a) sigma_K^{f*f}(t) via direct formula
  (b) sigma_K^{f*f}(t) reconstructed from the Fejer-Riesz coefficients
      used in v51's encoding
and verify they agree pointwise.
"""

import numpy as np


def sigma_direct(K: int, t_grid: np.ndarray) -> np.ndarray:
    """Direct computation of sigma_K^{f*f}(t) for uniform f = 2 * 1_{[-1/4, 1/4]}.

    For uniform f: a_k = 2 sin(pi k / 2) / (pi k), b_k = 0.
    So (f*f)^hat(k) = (a_k - i b_k)^2 = a_k^2 (real, since b_k = 0).
    sigma_K^{f*f}(t) = 1 + 2 sum_{k=1}^K w_k a_k^2 cos(2 pi k t)
    """
    vals = np.ones_like(t_grid)
    for k in range(1, K + 1):
        w_k = 1.0 - k / (K + 1)
        a_k = 2.0 * np.sin(np.pi * k / 2) / (np.pi * k)
        vals += 2.0 * w_k * (a_k**2) * np.cos(2.0 * np.pi * k * t_grid)
    return vals


def sigma_via_fr_coefs(K: int, t_grid: np.ndarray) -> np.ndarray:
    """Reconstruct sigma_K^{f*f}(t) from the coefficients used in v51's
    Fejer-Riesz encoding:
        Re(p_hat(k)) = -w_k * (M[k,k] - M[K+k, K+k]) = -w_k * (a_k^2 - b_k^2)
        Im(p_hat(k)) =  2 w_k * M[k, K+k]            = 2 w_k a_k b_k
    where p(t) = Omega - sigma(t). So:
        sigma(t) = Omega - p(t)
                 = Omega - [(Omega - 1) + 2 Re sum_k p_hat(k) e^{2 pi i k t}]
                 = 1 - 2 Re sum_k p_hat(k) e^{2 pi i k t}

    For uniform f (b_k = 0): p_hat(k) = -w_k a_k^2, so
        sigma(t) = 1 + 2 sum_k w_k a_k^2 cos(2 pi k t).
    Should match sigma_direct.
    """
    vals = np.ones_like(t_grid)
    for k in range(1, K + 1):
        w_k = 1.0 - k / (K + 1)
        a_k = 2.0 * np.sin(np.pi * k / 2) / (np.pi * k)
        b_k = 0.0
        # p_hat(k) = -w_k (a_k - i b_k)^2 = -w_k (a_k^2 - b_k^2 - 2i a_k b_k)
        p_hat_re = -w_k * (a_k**2 - b_k**2)
        p_hat_im = 2.0 * w_k * a_k * b_k  # = 0 for uniform
        # sigma(t) = 1 - 2 Re sum p_hat(k) e^{2 pi i k t}
        # Re(p_hat(k) e^{2 pi i k t}) = p_hat_re cos - p_hat_im sin
        vals -= 2.0 * (
            p_hat_re * np.cos(2.0 * np.pi * k * t_grid)
            - p_hat_im * np.sin(2.0 * np.pi * k * t_grid)
        )
    return vals


def f_star_f_uniform(t: np.ndarray) -> np.ndarray:
    """Exact (f*f)(t) for f = 2 * 1_{[-1/4, 1/4]}, t in [-1/2, 1/2]."""
    return np.where(np.abs(t) <= 0.5, 4.0 * np.maximum(0.5 - np.abs(t), 0.0), 0.0)


if __name__ == "__main__":
    K = 16
    t_grid = np.linspace(-0.5, 0.5, 1001)

    s_direct = sigma_direct(K, t_grid)
    s_via_fr = sigma_via_fr_coefs(K, t_grid)
    diff = np.max(np.abs(s_direct - s_via_fr))
    print(f"K={K}: max |sigma_direct - sigma_via_fr| = {diff:.2e}")

    # Also compare with exact (f*f)(t):
    exact = f_star_f_uniform(t_grid)
    print(f"K={K}: max sigma_direct = {s_direct.max():.6f}")
    print(f"       max (f*f) exact = {exact.max():.6f}  (should be 2.0)")
    print(f"       max |sigma_direct - (f*f)| = {np.abs(s_direct - exact).max():.4f}")

    # At what t does sigma peak?
    t_peak_direct = t_grid[np.argmax(s_direct)]
    t_peak_exact = t_grid[np.argmax(exact)]
    print(f"       peak of sigma at t = {t_peak_direct:.4f}")
    print(f"       peak of (f*f) at t = {t_peak_exact:.4f}")

    print("\nFor larger K:")
    for K in [32, 64, 96]:
        s_direct = sigma_direct(K, t_grid)
        s_via_fr = sigma_via_fr_coefs(K, t_grid)
        diff = np.max(np.abs(s_direct - s_via_fr))
        print(f"K={K}: max |sigma_direct - sigma_via_fr| = {diff:.2e}, max sigma = {s_direct.max():.6f}")
