# v53 — Additional localizer $\cos^3(2\pi x)$

## Goal

v51 uses one Putinar localizer $h_1(x) = \cos(2\pi x)$, which is
$\ge 0$ on $[-1/4, 1/4]$ and $\le 0$ on $[1/4, 3/4]$.  The
Schmüdgen/Putinar hierarchy also admits **powers** of this
localizer, since $\cos^{2k+1}(2\pi x)$ preserves the sign of
$\cos(2\pi x)$ and is therefore a valid support-specific localizer.

In v53 we add the **cubic** localizer
$$
h_3(x) \;:=\; \cos^3(2\pi x) \;=\; \tfrac{3}{4}\cos(2\pi x) + \tfrac{1}{4}\cos(6\pi x),
$$
which is a real trig polynomial of degree 3, nonneg on
$[-1/4, 1/4]$, negative on $(1/4, 3/4)$.

The measure $d\nu_3 := h_3(x)\,f(x)\,dx$ is a positive measure (for
$f \ge 0$ supported in $[-1/4, 1/4]$), so by Bochner, its
Hermitian Toeplitz $T_{\nu_3}$ is PSD.

## Fourier moments of $\nu_3$

Using the standard product-to-sum Fourier identity:
$$
\hat{\nu_3}(k)
\;=\; \int f(x)\,\bigl(\tfrac{3}{4}\cos(2\pi x) + \tfrac{1}{4}\cos(6\pi x)\bigr)\,e^{-2\pi i k x}\,dx
$$
$$
\;=\; \tfrac{3}{8}\bigl(\hat f(k-1) + \hat f(k+1)\bigr)
    + \tfrac{1}{8}\bigl(\hat f(k-3) + \hat f(k+3)\bigr).
$$

The corresponding real and imaginary parts:
$$
\mathrm{Re}\,\hat{\nu_3}(k)
= \tfrac{3}{8}(a_{|k-1|} + a_{|k+1|}) + \tfrac{1}{8}(a_{|k-3|} + a_{|k+3|}),
$$
$$
\mathrm{Im}\,\hat{\nu_3}(k)
= -\tfrac{3}{8}\bigl(\mathrm{sgn}(k-1) b_{|k-1|} + \mathrm{sgn}(k+1) b_{|k+1|}\bigr)
  - \tfrac{1}{8}\bigl(\mathrm{sgn}(k-3) b_{|k-3|} + \mathrm{sgn}(k+3) b_{|k+3|}\bigr).
$$
For the Toeplitz of size $K+1$, this needs Fourier moments of $f$
up to index $K+3$.

## Validity

Identical to v51, but with $h_3$ instead of $h_1$:
1. $\cos^3(2\pi x) \ge 0$ on $\mathrm{supp}\,f \subset [-1/4, 1/4]$.
2. Hence $\nu_3 = h_3 \cdot f \ge 0$ as a measure.
3. Bochner's theorem gives $T_{\nu_3} \succeq 0$ at any size.
4. At rank-1 feasible $(a, b)$, $T_{\nu_3}$ is computed exactly,
   and is PSD.  No rank-1 point is excluded.  $\square$

## Why this is not redundant with v51

v51 uses $h_1 = \cos$.  v53 uses $h_3 = \cos^3$.  These are
**different** localizers (not proportional), so they give
independent constraints.  The Fourier couplings are different:
v51 couples frequencies $k \pm 1$, v53 couples $k \pm 1$ and $k \pm 3$.

## Construction

v53 = v51 + additional localized Bochner $T_{\nu_3} \succeq 0$ of
size $K+1$, using Fourier moments of $f$ up to index $K+3$.

## Status

Awaiting theorist verdict.  The new constraint is a
textbook-Schmüdgen-style higher-power localizer in the
trigonometric Putinar hierarchy.
