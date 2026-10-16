# v46 — Extended Bochner–Toeplitz of degree 2K

## Goal

In v45 the Bochner–Toeplitz constraint is enforced at degree $K$
(matching the degree of the Fejér–Riesz $\sigma_K$).  But the actual
positive measure $f$ has Fourier coefficients at *all* frequencies,
and the Bochner–Toeplitz at higher degree is a strictly tighter
constraint on the same $(a_1, \dots, a_K)$ (since it requires the
existence of consistent $a_{K+1}, \dots$ that complete the
positive-measure structure).

In v46 we add *auxiliary* variables $a_{K+1}, \dots, a_{2K}$ and
$b_{K+1}, \dots, b_{2K}$, with their own cell bounds from $p$, and
enforce a $(2K+1) \times (2K+1)$ Bochner–Toeplitz constraint on
the full sequence $(a_0, a_1, \dots, a_{2K})$ minus
$i\,(b_0, b_1, \dots, b_{2K})$.

This **does not change the Fejér–Riesz constraint**, which still
acts on $\sigma_K^{f*f}$ (degree $K$) — the new auxiliaries only
serve to *tighten* the bound on $a_1..K, b_1..K$ via Bochner.

## Validity

The aux variables must satisfy:
1. **Cell bounds:** $a_k \in [\sum_j p_j \cos_\min^{(k)}, \sum_j p_j \cos_\max^{(k)}]$
   for all $k \in \{1, \dots, 2K\}$, with the same cell partition.
   These are valid (cell-wise rectangular bounds on $\int f \cos$).

2. **Bochner–Toeplitz $(2K+1) \times (2K+1)$:** the Hermitian
   Toeplitz with first row $(1, c_1, \dots, c_{2K})$, where
   $c_k = a_k - i b_k$, is PSD. Valid by Bochner's theorem applied
   to the (compactly supported) positive measure $f\,dx$.

3. **No coupling change:** the rest of the SDP is unchanged.
   The Fejér–Riesz constraint still uses only $a_1..K, b_1..K$.

The relaxation contains all rank-1 feasible $f$:
- A real $f \ge 0$ has well-defined Fourier coefs $a_k = \int f \cos$
  and $b_k = \int f \sin$ for *all* $k$. Setting the auxiliaries to
  these values yields a feasible point (cell bounds satisfied
  trivially, Bochner–Toeplitz at any size satisfied trivially).
- All other constraints are unchanged from v45 and still valid.

Hence v46 is a valid relaxation, and **strictly tighter** than v45
whenever the larger Bochner–Toeplitz cuts off some $(a_1..K, b_1..K)$
that fit a $(K+1)$-dimensional Bochner cone but not a
$(2K+1)$-dimensional one.

## Implementation

- Variables $a, b$ become length $2K+1$ instead of $K+1$.
- The moment matrix $M$ stays at size $2K+1$ (only on $a_1..K, b_1..K$),
  since Fejér–Riesz needs only those.
- Cell-bound coefficients computed for $k = 1, \dots, 2K$.
- Bochner–Toeplitz block becomes $2(2K+1) \times 2(2K+1)$ real PSD.

The dominant new cost is the $(2K+1)\times(2K+1)$ Bochner-Toeplitz PSD
(roughly $4 \times$ the v45 Bochner block) and additional cell bounds
($\sim 2 \times$).  Solver cost should scale by $\sim 4 \times$
relative to v45.

## Expected gain

Hard to predict without running.  At small $K$, the gain may be
negligible (Bochner is often nearly tight already).  At larger $K$,
the gap from rank-1 grows and the higher-degree Bochner may provide
$10^{-3}$–$10^{-4}$ tightening.

## Status

Awaiting theorist verdict.  No new mathematics — extension of v45's
Bochner argument from degree $K$ to degree $2K$.
