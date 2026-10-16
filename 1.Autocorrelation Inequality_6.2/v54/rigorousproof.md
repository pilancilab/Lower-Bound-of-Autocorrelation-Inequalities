# v54 — Additional localizer $\cos(2\pi x) - \cos^2(2\pi x)$

## Goal

Add a second, algebraically independent Putinar localizer
$$
h_2(x) \;:=\; \cos(2\pi x) - \cos^2(2\pi x) = \cos(2\pi x)(1 - \cos(2\pi x))
$$
to the v51 construction.

## Why $h_2$ is support-specific

Let $c = \cos(2\pi x)$.  Then $h_2 = c(1 - c)$.

- On $[-1/4, 1/4]$: $c \in [0, 1]$, so $c(1-c) \in [0, 1/4]$ — **nonneg**.
- On $(1/4, 3/4)$: $c \in [-1, 0)$, so $c(1-c) < 0$ (since
  $c < 0$ and $1-c > 0$) — **negative**.

So $h_2$ is nonneg exactly on $\{x : \cos(2\pi x) \ge 0\} = [-1/4, 1/4]$
(mod $\mathbb Z$), identical support as v51's $h_1 = \cos(2\pi x)$,
but the function itself is **algebraically different**.

## Why $h_2$ is not a linear combination of nonneg Bochner–summands

At first sight, $h_2 = h_1 - c^2$ appears to be a difference, and
$c^2 = \cos^2(2\pi x)$ is nonneg everywhere.  One might conclude
that $\text{Bochner}(h_2 \cdot f)$ is the difference of two PSDs,
hence not automatically PSD.

More precisely:
$$
T_{\nu_2} \;=\; T_{\nu_1} - T_{\nu_{c^2}}
$$
where $\nu_{c^2} = c^2 f$.  The right-hand side is a **difference**
of PSD matrices, which is not PSD in general.  So the constraint
$T_{\nu_2} \succeq 0$ is genuinely new — it cannot be derived from
individual PSD constraints on $T_{\nu_1}$ and $T_{\nu_{c^2}}$.

## Entries of $\nu_2$

Fourier expansion:
$$
h_2(x) = \cos(2\pi x) - \tfrac{1}{2}\bigl(1 + \cos(4\pi x)\bigr)
= -\tfrac{1}{2} + \cos(2\pi x) - \tfrac{1}{2}\cos(4\pi x).
$$
So the Fourier coefficients of $h_2$ (on $\mathbb T$) are
$\hat{h_2}(0) = -1/2$, $\hat{h_2}(\pm 1) = 1/2$, $\hat{h_2}(\pm 2) = -1/4$,
and zero elsewhere.

Hence
$$
\hat{\nu_2}(k) \;=\; \sum_m \hat{h_2}(m)\,\hat f(k - m)
\;=\; -\tfrac{1}{2}\hat f(k)
     + \tfrac{1}{2}\bigl(\hat f(k-1) + \hat f(k+1)\bigr)
     - \tfrac{1}{4}\bigl(\hat f(k-2) + \hat f(k+2)\bigr).
$$

In terms of $(a_k, b_k)$:
$$
\mathrm{Re}\,\hat{\nu_2}(k) = -\tfrac{1}{2} a_{|k|}
  + \tfrac{1}{2}(a_{|k-1|} + a_{|k+1|})
  - \tfrac{1}{4}(a_{|k-2|} + a_{|k+2|}),
$$
$$
\mathrm{Im}\,\hat{\nu_2}(k) = +\tfrac{1}{2}\,\mathrm{sgn}(k) b_{|k|}
  - \tfrac{1}{2}\bigl(\mathrm{sgn}(k-1) b_{|k-1|} + \mathrm{sgn}(k+1) b_{|k+1|}\bigr)
  + \tfrac{1}{4}\bigl(\mathrm{sgn}(k-2) b_{|k-2|} + \mathrm{sgn}(k+2) b_{|k+2|}\bigr).
$$
(Note: the sign on the first term flips from $-\mathrm{sgn}$ to
$+\mathrm{sgn}$ because $\hat{h_2}(0) = -1/2$, and
$\mathrm{Im}(\hat f(k)) = -\mathrm{sgn}(k) b_{|k|}$.)

For the Toeplitz of size $K+1$, we need $f$-moments up to index
$K+2$.

## Validity

Identical to v51, with $h_2$ replacing $h_1$:
1. $h_2(x) \ge 0$ for $x \in \mathrm{supp}\,f \subset [-1/4, 1/4]$.
2. $\nu_2 = h_2 \cdot f$ is a positive measure on $\mathbb T$.
3. $T_{\nu_2} \succeq 0$ by Bochner's theorem.  $\square$

## Construction

v54 = v51 (which already has $T_{\nu_1} \succeq 0$ and plain Bochner)
plus the additional PSD constraint $T_{\nu_2} \succeq 0$, which
requires Fourier moments up to index $K+2$.

## Status

Awaiting theorist verdict.  The localizer $h_2$ is algebraically
distinct from $h_1$ and from the trivial $h_0 = 1$, so the new
constraint is not in the span of the existing ones.
