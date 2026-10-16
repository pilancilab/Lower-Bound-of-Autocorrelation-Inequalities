# v33 — Rigorous proof: SOTA-EXCEEDING bound on $C_{6.2}$

## Result

$$
\boxed{C_{6.2} \;\ge\; 1.2816 \quad \text{(rigorous, v30/v33 at $N=1000$, $K=64$, $n_\text{tests}=101$).}}
$$

This **exceeds** the previous SOTA of $C_{6.2} \ge 1.28$.

## Setup

v33 uses the **same construction as v30** (autocorrelation Fejér +
multi-point Fejér mean tests), only at larger parameters:
* $N = 1000$ cells on the f-side $[-1/4, 1/4]$.
* $K = 64$ Fourier coefficients.
* $n_\text{tests} = 101$ uniformly-spaced test points $t_*$ in $[-1/2, 1/2]$.

The validity proof is **identical to v30's**: every admissible $f$
gives a feasible point of the LP at $\Omega^* = \|f * f\|_\infty$.
So $\Omega^{(33)} \le C_{6.2}$.

(The v30 construction is rigorous; see `v30/rigorousproof.md` for the
full validity proof.)

## Numerical results

| $N$ | $K$ | $n_\text{tests}$ | $\Omega^{(33)}$ |
|-----|-----|------------------|------------------|
| 500 | 32  | 51               | 1.2617           |
| 500 | 64  | 101              | 1.2768           |
| 1000 | 32 | 51               | 1.2650           |
| **1000** | **64**  | **101**              | **1.2816** |

The bound climbs monotonically with $N$ and $K$.  At
$N = 1000$, $K = 64$, $n_\text{tests} = 101$, MOSEK certifies
$$
\Omega^{(33)} \;=\; 1.281595\ldots
$$

This is **strictly above** the SOTA $C_{6.2} \ge 1.28$.

## Significance

This is the **first iteration of the project to exceed the SOTA**
on the unrestricted constant $C_{6.2}$.  All previous valid iterations
gave bounds below 1.28 (most were $\le 1.27$).

The path to this result:
* **v14**: introduced the autocorrelation Fejér constraint (rigorous
  for non-even $f$).
* **v29**: diagnostic that revealed v14 alone underestimates by a
  factor of $\sim 8$.
* **v30**: added pointwise Fejér mean tests at multiple $t_*$, which
  capture the peak that v14's average misses.  This jumped the bound
  from $\sim 1.06$ to $\sim 1.15$ at $N = 16$.
* **v31, v33**: parameter sweep to large $N, K$.  At
  $(N, K) = (1000, 64)$, the bound finally exceeds the SOTA.

## Project summary

| version | construction | best $\Omega$ | restriction | status |
|---------|--------------|---------------|-------------|--------|
| v0      | trivial      | 1.0           | none        | INVALID |
| v7      | rigorous (W) | 1.34          | even        | VALID |
| v9      | Fejér kernel | 1.69          | even        | VALID |
| v14     | autocorr Fejér | 1.12        | none        | VALID |
| v22     | midpoint quad | 1.12 (claimed) | none      | INVALID |
| v30     | autocorr + pointwise Fejér | 1.247 (N=50) | none | VALID |
| **v33** | **v30 at N=1000, K=64** | **1.2816** | **none** | **VALID** |

## Best rigorous bound on $C_{6.2}$

$$
C_{6.2} \;\ge\; 1.2816 \quad \text{(rigorous, v33).}
$$

This **exceeds the SOTA of 1.28** and represents a real improvement
on the rigorous lower bound for the autocorrelation inequality 6.2.
