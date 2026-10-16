# v6 — Rigorous proof of the asymmetry-slack constraints

## Setup

Same notation as v5: $f \ge 0$, $\operatorname{supp} f \subset [-1/4, 1/4]$,
$\int f = 1$, $F = f * f$, $L = 1/(4N)$, $p_j = \int_{I_j} f$,
$Q[j, k] = p_j p_k$.

v6 introduces a vector of slack variables $s_j \ge 0$ ($j = 1, \dots, 2N$)
intended to lift the asymmetry $(p_j - p_{2N+1-j})^{2}$.

## New constraints

### (S1) Per-index slack lower bound

$$
s_j \;\ge\; Q[j, j] + Q[2N+1-j, 2N+1-j] - 2\, Q[j, 2N+1-j].
\qquad (S1)
$$

**Proof.**  In the rank-one truth $Q[j,k] = p_j p_k$, the right-hand side
equals $p_j^{2} - 2 p_j p_{2N+1-j} + p_{2N+1-j}^{2} = (p_j - p_{2N+1-j})^{2} \ge 0$,
and we set $s_j := (p_j - p_{2N+1-j})^{2}$, so the inequality becomes
equality.  In particular it holds.  $\square$

### (S2) Sum upper bound on slacks

$$
\sum_{j=1}^{2N} s_j \;\le\; 2\,\operatorname{tr}(Q) - 2\sum_{j=1}^{2N} Q[j, 2N+1-j].
\qquad (S2)
$$

**Proof.**  In the rank-one truth, summing the equality version of (S1)
over $j$ gives
$$
\sum_j (p_j - p_{2N+1-j})^{2}
= 2\sum_j p_j^{2} - 2 \sum_j p_j p_{2N+1-j}
= 2\,\operatorname{tr}(Q) - 2\sum_j Q[j, 2N+1-j].
$$
With $s_j = (p_j - p_{2N+1-j})^{2}$ this is exactly the equality form of
(S2). $\square$

## Validity

For every admissible $f$, the assignment
$\Omega^* = \|f * f\|_\infty$, $p^* = (\int_{I_j} f)$,
$Q^*[j,k] = p_j^* p_k^*$, $s_j^* = (p_j^* - p_{2N+1-j}^*)^{2}$
is feasible: it satisfies (S1) and (S2) with equality and all v0/v1
constraints. $\square$

## Numerical observation

At $N = 8$, MOSEK still returns $\Omega^{(6)} = 1$.  The reason is that
the slack variables $s_j$ are FREE inside the relaxation: setting them
to large values absorbs the constraint without forcing $Q$ to deviate
from the diagonal.  Specifically, for the diagonal $Q = (1/(2N))\,I$
with uniform $p = (1/(2N))\,\mathbf 1$:

* (S1) requires $s_j \ge 1/(2N) + 1/(2N) - 0 = 1/N$.
* (S2) requires $\sum s_j \le 2 \cdot 1 - 2 \cdot 0 = 2$.

Both constraints are satisfied by choosing $s_j = 1/N$ (so $\sum s_j = 2$),
which is consistent with the rotated-SOC lower bound $s_j \ge 0$ (recall
$(p_j - p_{2N+1-j})^{2} = 0$ for uniform $p$).  Thus the relaxation has
"phantom" asymmetry slack that absorbs the would-be contradiction.

The fundamental obstruction is that the moment cone of admissible $f$'s
is the same as the convex hull of single-$f$ moments, which includes
mixtures of basis-vector measures.  Any convex relaxation that uses only
$(p, Q, \text{linear or PSD lifts})$ on the $f$-side cannot strictly
improve on $\Omega = 1$.

A genuinely tighter bound requires either
(i) restricting to even (or otherwise symmetric) $f$, as in v5; or
(ii) working entirely on the $F$-side with constraints that exclude
uniform $F$, which is the approach v7 will explore via Fourier moments
of $F$ tied to $f$ via the squared identity.
