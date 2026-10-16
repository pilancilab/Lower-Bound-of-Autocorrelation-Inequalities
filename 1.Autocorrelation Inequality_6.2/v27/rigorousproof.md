# v27 — Iteration note: variance constraint analysis

## Setup

v26's diagnostic showed that v14's LP-optimal $p$ is non-uniform
(edge-emphasized, symmetric).  v27 attempts to add an f-side
variance bound to restrict the LP from this regime.

## Attempted constraint

We tried adding $\sum_j p_j x_j^{2} \le 1/16$, the maximum variance
of an admissible $f$ (achieved at $f$ concentrated at $\pm 1/4$).

**Observation.**  This constraint is **automatic** for any admissible
$f$ (since $f$ supported on $[-1/4, 1/4]$ and any nonneg measure has
variance bounded by $(1/4)^{2}$ around its mean).  Adding it does
not restrict the feasible set further.

To get a tightening, the constraint would need to be of the form
$\sum_j p_j x_j^{2} \le V$ for $V < 1/16$, but such a $V$ is *not*
universally valid (the f at $\pm 1/4$ achieves $V = 1/16$).

## Conclusion

The variance constraint approach doesn't help.  The fundamental
issue identified in v26 — that v14's relaxation admits non-uniform
LP-optimal $p$ — cannot be resolved by adding more universal linear
constraints on $p$.

The path forward must use a *different* relaxation framework, e.g.,
one that penalizes non-uniformity in $p$ in a way consistent with
the convex relaxation principle.  This is left for future work.

## Validity

v27 introduces no actual new constraint (only a tautology).
$\Omega^{(27)} = \Omega^{(14)}$ at the same parameters.
