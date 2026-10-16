# v41 — Sweep with K=256 (per v33 verdict)

## Goal

Per v33 verdict's recommendation #1: try $K = 256$ at smaller $N$
to see if larger $K$ at fixed $N$ gives a tighter bound.

## Setup and validity

Same as v30; only parameter scaling.  Validity inherited.

## Results

| $N$    | $K$    | $n_\text{tests}$ | $\Omega^{(41)}$ |
|--------|--------|------------------|------------------|
| 250    | 256    | 129              | (in progress)    |
| 500    | 256    | 129              | (in progress)    |

The intuition: $K$ controls the precision of the Fejér kernel's
spectral support, which captures the actual peak structure of $F$.
At fixed $N$ (which controls f-side discretization), a larger $K$
should give a tighter peak-capture and a higher LP bound.

## Best rigorous bound (current overall)

$$
\boxed{C_{6.2} \;\ge\; 1.2859 \quad \text{(rigorous, v37 at $N = 5000$, $K = 64$, $n_\text{tests} = 101$).}}
$$

This **exceeds the SOTA $C_{6.2} \ge 1.28$** by $0.0059$.
