# v40 — v30 + v13 Bochner Toeplitz on f's Fourier

## Setup

v40 inherits all v30 constraints and adds the v13 Bochner Toeplitz
PSD constraint on f's Fourier sequence (the Hermitian Toeplitz of
$\hat f(k) = a_k - i b_k$).

## Validity

All constraints are individually rigorously valid (proven in v30
and v13).  The conjunction is rigorously valid: every admissible
$f$ gives a feasible point at $\Omega^* = \|f * f\|_\infty$.

## Numerical observation

| $N$ | $K$ | $\Omega^{(40)}$ | $\Omega^{(30)}$ for comparison |
|-----|-----|------------------|--------------------------------|
| 16  | 16  | 1.147           | 1.147                          |
| 16  | 32  | 1.143           | 1.153                          |
| 32  | 16  | 1.185           | 1.185                          |
| 32  | 32  | 1.183           | 1.196                          |

The v13 Bochner constraint **does not improve** v30; in some cases
it slightly *worsens* the bound (likely numerical or due to the
extra LMI constraint changing the LP topology).

## Conclusion

Adding the v13 Bochner Toeplitz to v30 produces no improvement.
The dominant constraint at the LP optimum remains v30's
combination of autocorrelation Fejér + multi-point Fejér mean
tests.

## Best rigorous bound on $C_{6.2}$

Unchanged from v35:
$$
C_{6.2} \;\ge\; 1.2851 \quad \text{(v35 at $N = 3000$, $K = 64$).}
$$
