# v18 — Rigorous proof: v14 + L³ constraint

## Setup

v18 builds on v14 (autocorrelation Fejér) and adds:

* F-side cell averages $W_m$ on intervals of $[-1/2, 1/2]$, linked
  to the (p, Q) lift via the v7 (W) "cells fully inside" inequality.
* The L³-norm constraint via Jensen and the chain
  $\int F^3 \le \Omega^2$.

## The L³ constraint chain

For $F \ge 0$, $\int F = 1$, $F \le \Omega$:
$$
F^3 \;=\; F \cdot F \cdot F
\;\le\; \Omega \cdot F \cdot F
\;\le\; \Omega^2 \cdot F.
$$
Integrating: $\int F^3 \le \Omega^2$.

By Jensen on the convex function $x \mapsto x^3$ (for $x \ge 0$),
applied to $F$ with cell-uniform measure:
$$
\frac{1}{L_F} \int_{B_m} F^3\,dt \;\ge\; \Bigl(\frac{1}{L_F}\int_{B_m} F\,dt\Bigr)^{\!3}
\;=\; W_m^{3}.
$$
So $\int_{B_m} F^3 \ge L_F W_m^{3}$ and summing over $m$:
$$
\int F^3 \;\ge\; L_F \sum_{m=1}^{4N} W_m^{3}.
$$

Combining,
$$
L_F \sum_{m=1}^{4N} W_m^{3} \;\le\; \int F^3 \;\le\; \Omega^2.
\tag{L3}
$$

## Lifting (L3) to a convex constraint

(L3) is "convex sum ≤ convex function of $\Omega$", which is not
directly convex.  We use auxiliary variables:

* $u_m \ge W_m^{3}$  (rotated power cone, convex).
* $\omega \ge \Omega^{2}$  (epigraph of square, convex).
* $L_F \sum u_m \le \omega$  (linear).

Each is convex; the conjunction is convex.

## Validity

For every admissible $f$, the rank-one truth gives feasible
$(\Omega^*, p, Q, W, u_m = W_m^{3}, \omega = (\Omega^*)^{2})$.  All
constraints are satisfied with equality at the rank-one truth.
Therefore $\Omega^{(18)} \le C_{6.2}$.

## Numerical observation

At small $(N, K)$, MOSEK returns the same value as v14 alone
(e.g. $\Omega^{(18)} = 1.033$ at $N = 8, K = 8$).  The L³ constraint
**does not bind** at the LP optimum because the LP picks $W$ close
to uniform (which minimizes $\sum W^{3}$).

For the L³ constraint to bind, we would need an additional constraint
that forces $W$ to be non-uniform — for example, a sharp lower
bound on $\max W = \Omega$ via a separate point evaluation.

The conclusion: v18's L³ constraint is mathematically valid and
correct, but does not give numerical improvement at small $(N, K)$.
The bottleneck is the LP's freedom to choose $W$ near-uniform.

## Conclusion

v18 is rigorously valid but produces the same numerical bound as
v14.  The v14 / v17 framework (autocorrelation Fejér + parameter
sweep) remains the dominant approach for the unrestricted constant.
