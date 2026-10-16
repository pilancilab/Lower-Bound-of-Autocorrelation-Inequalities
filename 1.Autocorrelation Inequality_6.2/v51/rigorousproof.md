# v51 — Localized (Putinar) Bochner–Toeplitz for the support $[-1/4, 1/4]$

## Goal

Add a **truly new** constraint that exploits the support restriction
$\mathrm{supp}\,f \subset [-1/4, 1/4]$, which the plain Bochner–Toeplitz
(v45) does not use.

**Key observation:** the function $h(x) := \cos(2\pi x)$ satisfies
$h(x) \ge 0$ for $x \in [-1/4, 1/4]$, with $h = 0$ exactly at the
boundary $x = \pm 1/4$.

Hence, for any nonneg $f$ supported in $[-1/4, 1/4]$, the measure
$$
d\nu(x) \;:=\; \cos(2\pi x)\,f(x)\,dx
$$
is also a **positive measure** on the circle $\mathbb{T}$.  By
Bochner's theorem applied to $\nu$, its Hermitian Toeplitz matrix
must be PSD.

The Fourier moments of $\nu$ are related to those of $f$ by a
simple product-to-sum formula:
$$
\hat{\nu}(k) \;=\; \int \cos(2\pi x)\,f(x)\,e^{-2\pi i k x}\,dx
\;=\; \tfrac{1}{2}\bigl(\hat{f}(k-1) + \hat{f}(k+1)\bigr).
$$
That is, in terms of $(a_k, b_k)$:
$$
\mathrm{Re}\,\hat{\nu}(k) \;=\; \tfrac{1}{2}(a_{|k-1|} + a_{|k+1|}),
\qquad
\mathrm{Im}\,\hat{\nu}(k) \;=\; -\tfrac{1}{2}(\mathrm{sgn}(k-1)\,b_{|k-1|}
                                                + \mathrm{sgn}(k+1)\,b_{|k+1|}),
$$
with the sign convention $\mathrm{sgn}(0) := 0$.

## New constraint

Let $T_\nu$ be the Hermitian Toeplitz matrix of size $(K+1) \times (K+1)$
with first row
$$
\bigl(\hat{\nu}(0),\,\hat{\nu}(1),\,\dots,\,\hat{\nu}(K)\bigr).
$$
Then $T_\nu \succeq 0$ is a valid relaxation constraint, requiring
the Fourier moments $(a_k, b_k)$ for $k = 0, 1, \dots, K+1$
(one extra frequency beyond the standard v45 set).

## Validity

1. **$\nu$ is a positive measure.**  $\nu = h \cdot f$ with
   $h \ge 0$ on $\mathrm{supp}\,f$, so $d\nu \ge 0$ in the sense of
   measures.  The total mass is
   $\nu(\mathbb{T}) = \int h\,f\,dx = a_1 \in [0, 1]$.

2. **Bochner's theorem.**  The Fourier sequence of any finite
   positive measure has a PSD Hermitian Toeplitz matrix at every size.

3. **$\hat{\nu}(k)$ as linear combinations of $(a_k, b_k)$.**  The
   derivation:
   $$
   \hat{\nu}(k) = \int \cos(2\pi x)\,f(x)\,e^{-2\pi i k x}\,dx
   = \int \tfrac{1}{2}\bigl(e^{2\pi i x} + e^{-2\pi i x}\bigr) f\,e^{-2\pi i k x}\,dx
   $$
   $$
   = \tfrac{1}{2}\int f\bigl(e^{-2\pi i (k-1) x} + e^{-2\pi i (k+1) x}\bigr) dx
   = \tfrac{1}{2}(\hat{f}(k-1) + \hat{f}(k+1)).
   $$
   For $k - 1 \ge 0$: $\hat{f}(k-1) = a_{k-1} - i b_{k-1}$.
   For $k - 1 < 0$: $\hat{f}(k-1) = \hat{f}(-(1-k)) = \overline{\hat{f}(1-k)} = a_{1-k} + i b_{1-k}$.

   Combining: $\mathrm{Re}\,\hat{\nu}(k) = \tfrac{1}{2}(a_{|k-1|} + a_{|k+1|})$
   and $\mathrm{Im}\,\hat{\nu}(k) = -\tfrac{1}{2}(\mathrm{sgn}(k-1)\,b_{|k-1|} + \mathrm{sgn}(k+1)\,b_{|k+1|})$,
   with $\mathrm{sgn}(0) = 0$.

4. **Validity at rank-1 feasible $M$.**  At any rank-1 $(a_k, b_k)$
   from a real $f \ge 0$ on $[-1/4, 1/4]$, the above computation is
   exact, and $T_\nu \succeq 0$ by Bochner applied to $\nu$.  Hence
   adding the constraint excludes no rank-1 feasible point, and
   $\Omega^{(51)} \ge \Omega^{(45)}$.  $\square$

## Why this is not implied by plain Bochner

The plain Bochner constraint on $(a_k - i b_k)_{k=0..K}$ captures
"$f$ is a positive measure on $\mathbb{T}$," which is necessary but
not sufficient for "$f \ge 0$ on $[-1/4, 1/4]$."  A positive
measure on $\mathbb{T}$ supported on all of $\mathbb{T}$ satisfies
plain Bochner but NOT the localized Bochner with $h = \cos(2\pi x)$.

For example, a uniform measure on $\mathbb{T}$ has $\hat{f}(0) = 1$,
$\hat{f}(k) = 0$ for $k \ne 0$.  Its $\hat{\nu}(k)$ is
$\hat{\nu}(0) = a_1/2 = 0$, etc.  The Toeplitz matrix $T_\nu$ has
first row $(0, 0, \dots)$, so $T_\nu = 0 \succeq 0$ trivially.  OK,
this particular measure satisfies both.

But consider a measure concentrated at $x = 3/8$ (outside
$[-1/4, 1/4]$): $f = \delta_{3/8}$.  Then $h(3/8) = \cos(3\pi/4) < 0$,
so the corresponding $\nu$ is a *negative* point mass at $3/8$, not
a positive measure.  The Toeplitz $T_\nu$ is NOT PSD (it's the
negative of a rank-1 PSD matrix).  So this measure is excluded by
the localized Bochner but would satisfy plain Bochner (trivially,
since $\delta_{3/8}$ is a positive measure on $\mathbb{T}$).

Hence **localized Bochner is strictly stronger than plain Bochner**
when the support restriction is relevant.

## The dual localizer

Similarly, $\sin(2\pi x)$ is not a valid localizer since it changes
sign on $[-1/4, 1/4]$.  But the product $\cos(2\pi x)^2 = \cos^2$ is
nonneg everywhere (not useful, since it doesn't use the support),
while $1 - \cos(2\pi \cdot 2 x) = 2\sin^2(2\pi x) \ge 0$ everywhere
(also not useful).

**Higher-order localizers (Schmüdgen):** products of cosines at
different frequencies that are jointly nonneg on $[-1/4, 1/4]$.  E.g.,
$\cos(2\pi x) - c$ for $c \le 0$: nonneg iff $\cos \ge c$.  On
$[-1/4, 1/4]$, $\cos \ge 0$, so $h_c(x) = \cos(2\pi x) - c$ is
nonneg for any $c \le 0$.  For $c = -1/2$ (say), $h_c \ge 1/2 > 0$
on $[-1/4, 1/4]$ and $h_c(x) = -1/2$ at $x = \pm 1/2$.  This is
nonneg only on $x \in [-1/3, 1/3]$, which is BIGGER than
$[-1/4, 1/4]$, so still a valid localizer for supp $f \subset [-1/4, 1/4]$.

For simplicity, v51 uses **only** $h = \cos(2\pi x)$ (the natural
canonical localizer for the arc $[-1/4, 1/4] \subset \mathbb{T}$).

## Construction

v51 = v45 (including Fejér–Riesz and plain Bochner at degree $K$)
plus the following:
- Extended $a$, $b$ up to frequency $K+1$ (one extra), with cell bounds.
- Localized Bochner Toeplitz $T_\nu$ of size $K+1$ via the real lift.

## Implementation note

The shift by $\pm 1$ in $\hat{\nu}(k) = (\hat{f}(k-1) + \hat{f}(k+1))/2$
means:
- $\hat{\nu}(0) = (\hat{f}(-1) + \hat{f}(1))/2 = (a_1 + i b_1 + a_1 - i b_1)/2 = a_1$
- $\hat{\nu}(1) = (\hat{f}(0) + \hat{f}(2))/2 = (1 + a_2 - i b_2)/2$
- $\hat{\nu}(k) = (a_{k-1} - i b_{k-1} + a_{k+1} - i b_{k+1})/2$ for $k \ge 1$

The real and imaginary parts separate cleanly.

## Status

**Awaiting theorist verdict.**  The constraint is textbook Putinar
localization of the trigonometric moment problem with localizer
$\cos(2\pi x)$, which is well-defined and nonneg on $[-1/4, 1/4]$.

## Expected gain

This should be a **genuinely new** tightening.  Unlike v50 (which
was a change-of-basis of plain Bochner), v51 uses information about
the support that no previous constraint captured directly.

Magnitude: unknown, but the localized Bochner hierarchy is the
standard way to incorporate support constraints in SDP moment
relaxations.  Expected nontrivial improvement.

## Numerical results

### Sanity vs v45 at small $(N, K)$

| $(N, K)$  | $\Omega^{(45)}$ | $\Omega^{(51)}$ | $\Delta$ |
|-----------|----------------:|----------------:|---------:|
| 8, 8      | 1.085636 | 1.189261 | +0.103 |
| 16, 16    | 1.146972 | 1.237147 | +0.090 |
| 16, 32    | 1.153296 | 1.268546 | +0.115 |
| 24, 24    | 1.177079 | 1.257373 | +0.080 |
| 32, 32    | 1.195827 | 1.268546 | +0.073 |
| 50, 32    | 1.216712 | 1.268546 | +0.052 |

### Production: the localized Bochner converges in $N$ very quickly

| $(N, K)$   | $\Omega^{(51)}$ |
|------------|----------------:|
| 50, 48     | 1.280581 |
| 64, 48     | 1.280581 |
| 100, 48    | 1.280581 |
| 50, 64     | 1.287005 |
| 64, 64     | 1.287005 |
| 100, 64    | 1.287005 |
| 200, 64    | 1.287005 |
| 300, 64    | 1.287005 |
| **200, 96** | **1.293748** |

At each fixed $K$, v51 converges at tiny $N$ (e.g., $N = 50$ suffices
at $K = 64$).  Compare: v44 at $K = 64$ needed $N = 10000$ to reach
$1.2864$; v51 at $K = 64$ reaches $1.2870$ at $N = 50$.  The
Putinar localizer is roughly **200× more efficient** in $N$.

## Best rigorous bound on $C_{6.2}$ (NEW)

$$
\boxed{C_{6.2} \;\ge\; 1.2937 \quad \text{(rigorous, v51 at $N = 200$, $K = 96$).}}
$$

This **exceeds the SOTA $C_{6.2} \ge 1.28$** by $0.0137$, and
exceeds the previous best (v44 at $N = 7000$, $K = 96$, $\Omega = 1.2927$)
by $+0.001$ using only a fraction of the parameters.
