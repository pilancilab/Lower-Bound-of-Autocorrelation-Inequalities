# v19 — Rigorous proof: v14 with multiple K Fejér constraints

## Setup

v19 inherits the v14 framework (autocorrelation Fejér via the SDP
lift on $(1, a, b)$) and adds the **same** Fejér-style constraint
at multiple values of $K$ simultaneously:
$$
\Omega \;\ge\; 1 + 2\sum_{k=1}^{K_\text{eff}}\Bigl(1 - \frac{k}{K_\text{eff}+1}\Bigr) v_k,
\qquad K_\text{eff} \in \{2, 4, 6, \dots, K_\text{max}\}.
$$

Each individual constraint is rigorously valid (proven in v14), and
adding several at different $K_\text{eff}$ only restricts the
feasible set further.

## Validity

Same as v14, applied at each $K_\text{eff}$.

## Numerical results

| $N$ | $K_\text{max}$ | $\Omega^{(19)}$ |
|-----|----------------|-----------------|
| 8   | 16             | 1.036           |
| 8   | 32             | 1.038           |
| 16  | 16             | 1.063           |
| 16  | 32             | 1.066           |
| 32  | 16             | 1.083           |
| 32  | 32             | 1.088           |

The improvement over v14 (which uses $K = K_\text{max}$ alone) is
small — about 0.001-0.003 — because the smaller-$K$ constraints
are dominated by the larger-$K$ one (Fejér weights at smaller $K$
include all the larger-$K$ frequencies with equal or greater
positive coefficients in many cases).

The dominant constraint remains the largest-$K$ Fejér.  v19's
multi-$K$ approach is rigorous but does not significantly improve.

## Conclusion

The autocorrelation Fejér framework (v14) gives bounds that increase
slowly with $N, K$ but appears to converge to the asymptotic limit
$\inf_f \int F^{2}$ which is at most $4/3$ (achieved at uniform $f$).
The convergence is slow, suggesting that genuinely new ideas are
needed to break past $\sim 1.10$ at moderate parameters.
