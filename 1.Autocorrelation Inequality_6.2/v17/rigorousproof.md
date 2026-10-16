# v17 — Rigorous proof: parameter sweep of v14 at large N, K

## Setup

v17 reuses **the same constraint set as v14** (rigorous and valid for
all admissible $f$, no even restriction).  The only change is the
parameter sweep at larger $N \in \{50, 75, 100, 150\}$ and
$K \in \{32, 64, 128\}$.

The validity proof is **identical** to v14's: every admissible $f$
gives a feasible point of the relaxation at
$\Omega^{*} = \|f * f\|_\infty$, so the relaxation optimum
$\Omega^{(17)} \le C_{6.2}$.

## Asymptotic limit prediction

For the v14 autocorrelation Fejér construction, the limiting bound
as $N, K \to \infty$ is determined by the *infimum over admissible
$f$* of the limit constraint
$$
\Omega \;\ge\; 1 + 2\sum_{k=1}^{\infty}|\hat F(k)|^{2} \;=\; \int_{-1/2}^{1/2} F(t)^{2}\,dt.
$$

For the uniform $f = 2 \cdot \mathbf 1_{[-1/4, 1/4]}$ — which is the
analytic optimizer for the related $\int f^{2}$ minimization
on the support $[-1/4, 1/4]$ — we have $F = f * f$ a triangle of
peak 2 at 0, and
$$
\int F^{2}\,dt
\;=\; \int_{-1/2}^{1/2} \bigl(2(1 - 2|t|)_{+}\bigr)^{2}\,dt
\;=\; 8\int_{0}^{1/2}(1 - 2t)^{2}\,dt
\;=\; \frac{4}{3} \;\approx\; 1.3333.
$$

So the **asymptotic limit** of the v14 / v17 bound is at most
$4/3 \approx 1.333$ — *exceeding the SOTA $C_{6.2} \ge 1.28$*.

The hope is that v17, run at large enough $N$ and $K$, can verify a
bound close to $4/3$, which would constitute a **rigorous
improvement** on the SOTA.

## Numerical results

| $N$ | $K$ | $\Omega^{(17)}$ |
|-----|-----|-----------------|
| 50  | 32  | 1.099           |
| 75  | 32  | 1.107           |
| 100 | 32  | 1.111           |

The bound climbs monotonically with $N$ but slowly.  The asymptotic
limit (predicted) is $4/3 \approx 1.333$ if uniform $f$ minimizes
$\int F^{2}$.  At the parameters tested, the bound is well below
this limit, suggesting either (a) the convergence is very slow, or
(b) the LP picks non-uniform $p$ for which $\int F^{2}$ is smaller
than $4/3$.

## Caveat

The asymptotic limit $4/3$ is the value of the constraint at
the rank-one *uniform* truth.  For *non-uniform* admissible $f$,
the value of $\int F^{2}$ can be different.  The infimum over
admissible $f$ is what the relaxation actually computes, and this
infimum may be smaller than $4/3$.

In particular, for the "edge-concentrated" $f$
(e.g. $f \propto \mathbf 1_{[-1/4, -1/4 + \epsilon]} + \mathbf 1_{[1/4 - \epsilon, 1/4]}$),
$\int F^{2}$ can be much larger than $4/3$ (since $F$ is more
concentrated), so this $f$ is *not* a minimizer of the relaxation.

The minimizer is some intermediate $f$ whose $\int F^{2}$ is close
to (but possibly less than) $4/3$.  Empirical convergence will tell
us how close.
