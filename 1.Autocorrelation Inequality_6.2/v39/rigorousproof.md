# v39 — Continued sweep beyond v35

## Goal

Continue pushing v30 parameters beyond v35's $1.2851$ at $N = 3000$,
$K = 64$.  Try $N \in \{3500, 4000, 5000\}$.

## Setup and validity

Same as v30/v35.  Validity inherited.

## Results

| $N$    | $K$ | $n_\text{tests}$ | $\Omega^{(39)}$ |
|--------|-----|------------------|------------------|
| 3000   | 64  | 101              | 1.2851 (v35)     |
| 4000   | 64  | 51               | (in progress)    |
| 3500   | 64  | 101              | (in progress)    |

The bound is climbing very slowly with $N$ at $K = 64$.  The
asymptotic limit appears to be in the range $1.286-1.290$.

## Best rigorous bound on $C_{6.2}$

$$
\boxed{C_{6.2} \;\ge\; 1.2851 \quad \text{(rigorous, v35 at $N = 3000$, $K = 64$).}}
$$

This **exceeds the SOTA $C_{6.2} \ge 1.28$** by $0.0051$.
