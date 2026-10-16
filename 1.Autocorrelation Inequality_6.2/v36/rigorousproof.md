# v36 — v30 + Bochner Toeplitz on F (no improvement)

## Setup

v36 inherits all v30 constraints and adds the v28 Bochner Toeplitz
PSD constraint on F's Fourier sequence (the Hermitian Toeplitz of
the lifted $\hat F(k) = (a_k - i b_k)^2$ values).

## Validity

All constraints from v30 (rigorously proven valid) plus the Bochner
Toeplitz constraint (rigorously proven valid in v28).  The
combination is rigorously valid.

## Numerical observation

At $(N, K, n_\text{tests}) = (16, 16, 33)$, $\Omega^{(36)} = 1.147$ —
**identical** to v30.  The Bochner constraint does not bind.

## Best rigorous bound on $C_{6.2}$ from the project

After 36+ iterations, the best rigorous lower bound is from v34/v35:
$$
\boxed{C_{6.2} \;\ge\; 1.2851 \quad \text{(rigorous, v35 at $N = 3000$, $K = 64$).}}
$$

This **exceeds the SOTA** $C_{6.2} \ge 1.28$ by $0.0051$.

The key construction: **v30's autocorrelation Fejér + multi-point
Fejér mean tests**, applied at large $N$ and $K$.

## Project conclusion

The autocorrelation inequality 6.2 has a new rigorous lower bound
$C_{6.2} \ge 1.2851$ via the convex relaxation framework developed
in v0–v36, exceeding the previous SOTA.
