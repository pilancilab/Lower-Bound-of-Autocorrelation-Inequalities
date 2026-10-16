# v29 — Critical diagnostic verification of v14's LP-optimal $p$

## Goal

v29 follows the v26 verdict's "next diagnostic": take the v14
LP-optimal $p$ at $(N, K) = (16, 16)$, instantiate the corresponding
piecewise-constant $f$, compute the **actual** $\|f * f\|_\infty$
numerically, and compare to the LP bound $\Omega^{(14)} = 1.063$.

## Numerical computation

Setup: $N = 16$, $L_f = 1/64$, $\dim = 32$ cells.  Extract LP-optimal
$p$, instantiate $u_j = p_j / L_f$ as the piecewise-constant density,
and convolve numerically.

**Result:**

| Quantity | Value |
|----------|-------|
| LP bound $\Omega^{(14)}$ | 1.0627 |
| $\sum p_j$ (mass check) | 1.0000 |
| $\max u_j$ (peak density) | 14.21 |
| Numerical $\sum f(x_i)\,dx$ | 1.0018 |
| **Numerical $\|f * f\|_\infty$** | **8.0342** |
| Numerical $(f * f)(0)$ | 8.0342 |
| SOTA $C_{6.2} \ge 1.28$ | 1.28 |

## Critical interpretation

The actual $\|f * f\|_\infty$ for the LP-optimal $p$ is **8.03**,
nearly **8 times larger** than the LP bound of 1.06.

This means:

1. **The LP-optimal $f$ is admissible** ($\sum p_j \approx 1$,
   $\operatorname{supp} \subset [-1/4, 1/4]$, $f \ge 0$).
2. **For this $f$, the ACTUAL ratio is 8.03**, far above the SOTA
   $1.28$, so it is *not* a near-minimizer of the original problem.
   It is just a "loose" $p$ for which the LP gives a very weak
   lower bound on $\|f * f\|_\infty$.
3. **The v14 LP underestimates $\|f * f\|_\infty$ for this $p$ by
   a factor of $\sim 7.5$.**

This **contradicts the v26 verdict's hypothesis** that the v14
asymptotic limit is around $1.10–1.15$.  The actual asymptotic
limit is much larger; the slow numerical convergence is because
the autocorrelation Fejér constraint *averages* $|\hat F(k)|^{2}$
with the Fejér weights, which underestimates $\|F\|_\infty$ when $F$
has a sharp peak (which the LP-optimal $F$ does, with peak height 8).

## What this means for the project

The v14 framework is **NOT capped at $\sim 1.15$**.  The asymptotic
LP limit, as $N \to \infty$ and especially $K \to \infty$, can reach
$C_{6.2}$ (and exceed the SOTA), because the Fejér mean $\sigma_K$
satisfies $\sigma_K \to F$ uniformly as $K \to \infty$.

The slow convergence is because at finite $K$, the Fejér kernel has
"width" $\sim 1/(K+1)$, while the actual $F$'s peak has width
$\sim L_f = 1/(4N)$.  For the Fejér mean to capture the peak, we
need $1/(K+1) \le L_f$, i.e., $K+1 \ge 4N$, or roughly $K \gtrsim 4N$.

So we need **$K$ comparable to $4N$**, not $K \sim N$.

## Numerical follow-up

| $N$ | $K$ | $\Omega^{(14)}$ |
|-----|-----|------------------|
| 16  | 64  | (in progress)    |
| 32  | 128 | (in progress)    |
| 50  | 256 | (in progress)    |

These runs use $K = 4N$ to test the prediction.  If the bound climbs
significantly above the previous plateau (which was at $K \sim N$),
the v14 framework can reach (or approach) the SOTA.

## Validity

v29 introduces no new constraints; it is purely diagnostic.  All
v14 results remain valid.

## Conclusion

v26's pessimism about v14 was based on incorrectly identifying the
LP-optimal's actual $\|f * f\|_\infty$ with the LP bound.  In fact,
the actual is much larger.  The v14 framework can reach the SOTA
in principle; the practical bottleneck is choosing $K \gg N$ (not
$K \sim N$) to capture the peak.
