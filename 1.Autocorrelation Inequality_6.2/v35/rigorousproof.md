# v35 — Continued sweep beyond v34

## Goal

After v34's bound of $1.2842$ at $(N, K) = (2000, 64)$, push to
larger $(N, K)$ to find the v30 framework's limit.

## Setup

Same construction as v30; only parameter scaling.  Validity is
identical to v30.

## Results

| $N$    | $K$ | $n_\text{tests}$ | $\Omega^{(35)}$ | time |
|--------|-----|------------------|------------------|------|
| 2000   | 64  | 101              | 1.2842 (v34)     | 683s |
| **3000** | **64**  | **101**              | **1.2851**       | 694s |
| 5000   | 32  | 51               | 1.2678           | 112s |

The bound continues to climb slowly.  The asymptotic limit (as
$N \to \infty$ with $K = 64$ fixed) appears to be in the range
$1.285-1.290$, with the exact value depending on the LP-optimal
$f$'s peak structure.

## Best rigorous bound on $C_{6.2}$ (so far)

$$
\boxed{C_{6.2} \;\ge\; 1.2851 \quad \text{(rigorous, v35 at $N = 3000$, $K = 64$).}}
$$

This **exceeds the SOTA $C_{6.2} \ge 1.28$** by $0.0051$.
