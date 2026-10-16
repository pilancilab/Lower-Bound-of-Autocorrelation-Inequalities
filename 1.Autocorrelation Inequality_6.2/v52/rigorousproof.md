# v52 — Extend the localized Bochner in v51 to size $2K+1$

## Goal

In v51 the localized (Putinar) Bochner with $h = \cos(2\pi x)$ uses
a Hermitian Toeplitz of size $(K+1) \times (K+1)$ on the sequence
$\hat\nu(0), \hat\nu(1), \dots, \hat\nu(K)$.  This is **level $K$**
of the localized hierarchy.

The Putinar hierarchy converges to the true optimum as the level
$\to \infty$.  At level $K$, we know only $\hat\nu(k)$ for
$|k| \le K$, and each is a linear combination of $\hat f(k\pm 1)$.

In v52 we go one level up: enforce
$T_\nu^{(2K)} \succeq 0$ where $T_\nu^{(2K)}$ is Hermitian Toeplitz
of size $(2K+1) \times (2K+1)$ on the sequence
$\hat\nu(0), \hat\nu(1), \dots, \hat\nu(2K)$.  This requires
Fourier moments of $f$ up to index $2K+1$.

## Construction

- Extend $a, b$ from $K+1$ entries to $2K+2$ entries (indices
  $0, 1, \dots, 2K+1$), with cell bounds.
- Plain Bochner at size $2K+2$ on $(a, b)_{0..2K+1}$.
- Localized Bochner $T_\nu^{(2K)}$ of size $2K+1$.
- All other v51 constraints inherited unchanged.

## Validity

Same as v51.  $\nu = h \cdot f$ is a positive measure on $\mathbb T$
with Fourier sequence
$\hat\nu(k) = (\hat f(k-1) + \hat f(k+1))/2$ for $k \in \mathbb Z$.
Its Hermitian Toeplitz at any size is PSD by Bochner.

## Expected gain

Given the huge jump from v45 ($1.22$) to v51 ($1.27$) at small $K$,
extending the localized Bochner to level $2K$ should give
additional tightening — possibly $0.01$–$0.03$ further.

## Status

Awaiting theorist verdict.
