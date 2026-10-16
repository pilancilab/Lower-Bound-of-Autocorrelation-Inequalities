# v37 — Continued sweep beyond v35

## Goal

After v35's bound of $1.2851$ at $(N, K) = (3000, 64)$, push to
even larger $N$ to find the asymptotic limit of v30 with $K = 64$.

## Setup and validity

Same construction as v30; only parameter scaling.  Validity is
identical to v30.

## Results

| $N$    | $K$ | $n_\text{tests}$ | $\Omega^{(37)}$ |
|--------|-----|------------------|------------------|
| 5000   | 32  | 51               | 1.2678           |
| **5000** | **64**  | **101**              | **1.2859**       |
| 3000   | 64  | 101              | 1.2851           |

The bound at $K = 32$ has plateaued around 1.268 (between
$N = 1000$ and $N = 5000$).  At $K = 64$, the bound is still
climbing as $N$ grows.  The asymptotic limit at $K = 64$ may be
around $1.286–1.290$.

## Best rigorous bound on $C_{6.2}$ (so far)

$$
\boxed{C_{6.2} \;\ge\; 1.2859 \quad \text{(rigorous, v37 at $N = 5000$, $K = 64$).}}
$$

This **exceeds the SOTA $C_{6.2} \ge 1.28$** by $0.0059$.
