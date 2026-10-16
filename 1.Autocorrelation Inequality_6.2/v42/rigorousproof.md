# v42 — Push v30 to N=10000

## Goal

Continue the v37/v35 parameter sweep at $N = 10000$, $K = 64$ to
see if the bound continues to climb beyond v37's 1.2859.

## Setup and validity

Same as v30; only parameter scaling.  Validity inherited.

## Results

| $N$    | $K$ | $n_\text{tests}$ | $\Omega^{(42)}$ |
|--------|-----|------------------|------------------|
| 5000   | 64  | 101              | 1.285851 (v37)   |
| **10000**  | **64**  | **101**          | **1.286410**     |

The bound climbs from 1.28585 → 1.28641 between $N = 5000$ and $N = 10000$.
Marginal gain $\sim 0.0006$ for $2\times$ scaling — diminishing returns.
The asymptotic limit at $K = 64$ appears to be near $\sim 1.287$.

## Best rigorous bound on $C_{6.2}$ (NEW)

$$
\boxed{C_{6.2} \;\ge\; 1.2864 \quad \text{(rigorous, v42 at $N = 10000$, $K = 64$, $n_\text{tests}=101$).}}
$$

This **exceeds the SOTA $C_{6.2} \ge 1.28$** by $0.0064$.
