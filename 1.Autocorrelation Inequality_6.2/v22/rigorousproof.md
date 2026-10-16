# v22 — Proof attempt and self-INVALID acknowledgement

## Goal

v22 attempts to tighten the cell-bounds on $a_k$, $b_k$ in v14 by
replacing the loose min/max-of-cos cell bound with a midpoint
quadrature plus an explicit error.  The hope is to accelerate the
convergence of the LP bound to its asymptotic limit.

## Attempted constraint

For the midpoint approximation with error
$$
\Bigl|\,a_k - \sum_{j} p_j\cos(2\pi k\,x_j)\,\Bigr|\;\le\;\mathrm{err}_k,
$$
we tried $\mathrm{err}_k = \pi^{2}k^{2}L_f^{2}/2$, motivated by the
*standard midpoint rule* error for smooth integrands.

## Why this is INVALID

The standard midpoint rule error
$\mathrm{err}_k \sim L_f^{2}\cdot\sup|g''|$ assumes the integrand
is smooth (specifically, $g''$ exists and the *uniform measure* is
integrated against $g$).  In our case the integrand is $f \cdot g$
where $f$ is **only** $L^{1}$ (a nonneg density on $I_j$ with mass
$p_j$, NOT a smooth function).  The midpoint rule on a non-smooth
integrand has error $O(L_f)$, not $O(L_f^{2})$.

**Concrete counterexample.**  Take $f = \delta_{x^{*}}$ (concentrated
at a point $x^{*}$) on cell $I_j$.  Then $p_j = 1$ and
$$
\int_{I_j} f(x)\cos(2\pi k x)\,dx \;=\; \cos(2\pi k x^{*}),
$$
while the midpoint approximation gives $\cos(2\pi k\, x_j)$ where
$x_j$ is the cell midpoint.  The error is
$|\cos(2\pi k x^{*}) - \cos(2\pi k x_j)| \le |2\pi k|\cdot|x^{*} - x_j|
\le (2\pi k)\cdot(L_f/2) = \pi k L_f$,
which is $O(L_f)$, not $O(L_f^{2})$.

For $N = 8$, $k = 1$, $L_f = 1/32$: the *true* error can be as large
as $\pi/32 \approx 0.098$, but v22 only allows $\pi^{2}/2 \cdot (1/32)^{2}
\approx 0.0048$.  The constraint is *too tight* by a factor of about
$20$, and **excludes the concentrated $f$**, which is admissible.

So v22's relaxation may **not** contain every admissible $f$ as a
feasible point, and the LP optimum is **not** certifiably $\le C_{6.2}$.

## Numerical observation

| $N$ | $K$ | $\Omega^{(22)}$ |
|-----|-----|-----------------|
| 8   | 8   | 1.087           |
| 8   | 16  | 1.096           |
| 8   | 32  | 1.100           |
| 16  | 8   | 1.092           |
| 16  | 16  | 1.109           |
| 16  | 32  | 1.116           |
| 32  | 8   | 1.090           |
| 32  | 16  | 1.113           |
| 32  | 32  | 1.124           |

These are *higher* than v14 at the same parameters, which suggests
that the v22 constraint is genuinely tighter — but **the
tightening is non-rigorous**.  We cannot certify that
$\Omega^{(22)} \le C_{6.2}$.

## Verdict (self-INVALID)

v22 is **mathematically invalid** as written.  The error bound for
the midpoint quadrature on a non-smooth integrand is $O(L_f)$, not
$O(L_f^{2})$.  Replacing $\pi^{2}k^{2}L_f^{2}/2$ with the rigorous
$\pi k L_f$ recovers the loose v14 cell-bound (no tightening at all).

## v23 plan

Fix v22 by either:
1. **Using the rigorous $\pi k L_f$ error.**  This is just v14's loose
   cell-bound in disguise; no tightening.
2. **Adding a constraint that $f$ is bounded** ($\|f\|_\infty \le M$
   for some $M$ that the LP can compute), which would justify a
   higher-order quadrature error.  However, $\|f\|_\infty$ is *not
   bounded above* on the admissible set (delta-like $f$ have
   infinite $\|f\|_\infty$).

Option 2 doesn't work for arbitrary admissible $f$, but it could
work if we *restrict* to $f$ bounded by a known $M$ — at the cost
of bounding only $C_{6.2}^{M}$ (the bounded-$M$ version), which
$\le C_{6.2}$ in general.

v23 will explore this restricted bound.
