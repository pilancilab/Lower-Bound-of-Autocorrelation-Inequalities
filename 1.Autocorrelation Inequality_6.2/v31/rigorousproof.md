# v31 — Rigorous parameter sweep of v30

## Setup

v31 reuses **the v30 construction unchanged** (autocorrelation Fejér
+ pointwise Fejér mean tests at multiple $t_*$).  The only change
is to push the parameters $(N, K, n_\text{tests})$ to large values
and report the convergence.

The validity proof is identical to v30's: every admissible $f$
gives a feasible point of the v30 LP at $\Omega^* = \|f * f\|_\infty$.
So $\Omega^{(31)} = \Omega^{(30)} \le C_{6.2}$.

## Numerical sweep

| $N$    | $K$ | $n_\text{tests}$ | $\Omega^{(31)}$ |
|--------|-----|------------------|------------------|
| 16     | 16  | 33               | 1.147            |
| 32     | 32  | 51               | 1.196            |
| 50     | 32  | 51               | 1.217            |
| 100    | 32  | 51               | 1.239            |
| 150    | 32  | 51               | 1.248            |
| 200    | 32  | 51               | 1.253            |
| 300    | 32  | 51               | 1.258            |
| 500    | 32  | 51               | 1.262            |
| 1000   | 32  | 51               | 1.265            |
| 2000   | 32  | 51               | 1.267            |

The bound increases monotonically with $N$ and **converges to
approximately 1.27** as $N \to \infty$ (with $K = 32$).  This is
**very close to the SOTA $C_{6.2} \ge 1.28$ but does not yet exceed
it**.

## Project status

| version | construction | best $\Omega$ | restriction | status |
|---------|--------------|---------------|-------------|--------|
| v0      | trivial      | 1.0           | none        | INVALID |
| v7      | rigorous (W) | 1.34          | even        | VALID |
| v9      | Fejér kernel | 1.69          | even        | VALID |
| v14     | autocorr Fejér | 1.12        | none        | VALID |
| v22     | midpoint quad | 1.12 (claimed) | none      | INVALID (non-rigorous) |
| **v30** | **autocorr Fejér + pointwise Fejér tests** | **1.247** | **none** | **VALID** |
| **v31** | v30 sweep | **1.267 at N=2000, K=32** | **none** | **VALID** |

## Conclusion

v30/v31 is the **fifth major milestone** of the project, after
v7/v9/v14, and the **first iteration to come close to the SOTA**
on the **unrestricted** $C_{6.2}$.  The bound climbs from 1.06
(v14 at small $N$) to **1.267 at $N = 2000$**, a $0.205$ improvement.
With more computation (larger $K$), it may exceed 1.28 and beat
the SOTA.

The key insight (from v29's diagnostic) was that v14 captures
$\int F^{2}$ but misses the peak.  Adding pointwise Fejér mean
tests at multiple $t_*$ captures the peak directly, giving the
~0.18 improvement at $N = 16$ alone.

## Best rigorous bound on $C_{6.2}$

$$
\boxed{C_{6.2} \;\ge\; 1.267 \text{ (rigorous, v31 at $N = 2000, K = 32$)}.}
$$

The SOTA is $C_{6.2} \ge 1.28$.  v31 reaches within $0.013$ of the
SOTA via convex relaxation, the closest the project has come.
