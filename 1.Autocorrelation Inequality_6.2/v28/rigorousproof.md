# v28 — Bochner Toeplitz on F's Fourier sequence

## Setup

v28 inherits all v14 constraints and adds the **Bochner Toeplitz**
PSD constraint on the lifted F's Fourier sequence:
$$
T_{F}[i, j] := \hat F(j - i) \quad\text{(Hermitian Toeplitz, $i, j = 0, \dots, K$)},
$$
with $\hat F(0) = 1$ and $\hat F(k) = (a_{k}^{2} - b_{k}^{2}) - 2 i\, a_{k} b_{k}$
lifted via $M$:
$$
\operatorname{Re}\hat F(k) = M[k, k] - M[K{+}k, K{+}k],
\quad \operatorname{Im}\hat F(k) = -2\, M[k, K{+}k].
$$
The constraint $T_{F} \succeq 0$ (Hermitian PSD) is expressed as a real
block LMI.

## Validity

**Bochner's theorem.**  $F = f * f \ge 0$ implies the Hermitian
Toeplitz $[\hat F(j-i)]_{i,j}$ is positive semidefinite at every
size.  In particular at size $K + 1$.

**Validity in the lift.**  The rank-one $M = vv^{\top}$ from any
admissible $f$ gives the lifted $\hat F(k)$ matching the actual
$\hat F(k) = \hat f(k)^{2}$, so $T_{F}$ is the actual Toeplitz of
$F$'s Fourier coefficients, which is PSD by Bochner.  In the
relaxation (with $M$ not necessarily rank-one), the lifted $T_{F}$
might not be PSD; adding $T_{F} \succeq 0$ as a constraint just
restricts the feasible set further, preserving validity.  $\square$

## Numerical observation

| $N$ | $K$ | $\Omega^{(28)}$ |
|-----|-----|-----------------|
| 8   | 8   | 1.033           |
| 8   | 16  | 1.036           |
| 16  | 8   | 1.055           |
| 16  | 16  | 1.063           |

Identical to v14 at the same parameters.  The Bochner Toeplitz
constraint is **not binding** in v14's feasible region.

The reason is that v14's LP-optimal $M$ already has lifted $T_{F}$
PSD (because the LP uses the rank-one truth $M = vv^{\top}$ via
the $u_k = M[k,k] + M[K{+}k, K{+}k]$ construction, which forces
$T_{F}$ to be PSD automatically).

## Conclusion

v28 is **rigorously valid** but produces **no improvement** over v14.
The Bochner constraint on $F$'s Fourier sequence is implied by the
v14 SOC chain.

This is the **last constraint type** I attempted to add to v14
without sacrificing rigor.  All Fourier-side and PSD-side
restrictions have been tried; none improves the bound.

The path forward (per v20, v26 verdicts) is to either:
1. Run v14 at very large $N, K$ to approach its asymptotic limit
   (which v26 diagnostics suggest is around $1.10-1.15$, *not*
   the SOTA $1.28$).
2. Find a completely different relaxation framework (e.g., a
   sum-of-squares hierarchy on the variational problem of $f$).

Both are open work for future iterations.
