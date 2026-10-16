# v34 — Rigorous parameter sweep beyond v33

## Goal

v33 achieved the SOTA-exceeding bound $\Omega \ge 1.2816$ at
$(N, K, n_\text{tests}) = (1000, 64, 101)$.  v34 pushes the
parameters further to find the largest rigorous bound v30 can give.

## Setup

Same construction as v30/v31/v33; only parameter scaling.  Validity
is identical.

## Numerical results

| $N$    | $K$ | $n_\text{tests}$ | $\Omega^{(34)}$ | time |
|--------|-----|------------------|------------------|------|
| 1000   | 64  | 101              | 1.2816 (v33)     | 494s |
| 2000   | 64  | 101              | **1.2842**       | 683s |
| 1000   | 32  | 201              | 1.2651           | 136s |

**Best v34 bound: $C_{6.2} \ge 1.2842$ at $(N, K) = (2000, 64)$.**

This is **strictly above the SOTA $1.28$** by $0.0042$.

Increasing $n_\text{tests}$ from 51 to 201 at fixed $(N, K) = (1000, 32)$
gave essentially no improvement (1.265 vs 1.265), suggesting that
the test-point density is already saturated at $n_\text{tests} \approx 51$.
The binding parameters are $N$ and $K$.

The expected behavior is monotone increase, but the convergence rate
slows for $N \gg K$ (the f-side discretization is no longer the
bottleneck).

The key parameter to push appears to be **$K$ at fixed
moderate $N$**.  At $N \approx 1000$, $K = 64$ gives 1.2816;
$K = 96$ may push to $\sim 1.29$ or higher.

## Validity

Same as v30; rigorous lower bound on the unrestricted $C_{6.2}$.
