# v43 — Continue iterating with intermediate K values

## Setup

Same v30 construction; sweep at intermediate $K = 80, 96$ at moderate $N$.

## Validity

Inherited from v30.

## Results

| $N$    | $K$ | $n_\text{tests}$ | $\Omega^{(43)}$ |
|--------|-----|------------------|------------------|
| 500    | 96  | 101              | 1.2803           |
| 1000   | 96  | 101              | **1.2860**       |

The $(N=1000, K=96)$ run beats $(N=5000, K=64)$ by $\sim 0.00015$
and is achieved with $\sim$ half the f-side dimension.  This confirms
v33's recommendation: at this scale, **larger $K$ matters more than
larger $N$**.

## Best rigorous bound on $C_{6.2}$

Subsequent v42 result $(N=10000, K=64) \to 1.2864$ exceeds v43's
$(1000, 96) \to 1.2860$, so the overall best is unchanged from v42:
$$
C_{6.2} \;\ge\; 1.2864 \quad \text{(rigorous, v42 at $N = 10000$, $K = 64$).}
$$

The v43 result remains useful as a data point: it shows that the
$K = 96$ track at $N = 1000$ matches the $K = 64$ track at $N = 10000$
within $\sim 5 \cdot 10^{-4}$, suggesting both tracks converge to a
common limit $\sim 1.287$.
