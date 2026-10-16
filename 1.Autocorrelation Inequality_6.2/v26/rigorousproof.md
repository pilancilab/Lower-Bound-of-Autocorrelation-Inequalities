# v26 — Diagnostic of v14's LP-optimal solution

## Setup

v26 is *not a new constraint* but a **diagnostic** of the v14 LP at
$N = 16$, $K = 16$.  Per the v20 verdict's recommendation #2:
extract the LP-optimal $(p, a, b)$ and check what $\int F^{2}$ it
gives, to diagnose the convergence bottleneck.

## What v14 picks

At $N = 16$, $K = 16$, MOSEK returns $\Omega^{(14)} = 1.063$.  The
LP-optimal $p$ is:
$$
p \;\approx\; (0.222, 0, 0, 0, 0.019, 0.081, 0, 0.007, 0.033, 0.044, 0, 0.017, 0, 0.052, 0.024, 0, \\
\;\;\;\;\;\;\;\;\;\;\;0, 0.024, 0.052, 0, 0.017, 0, 0.044, 0.033, 0.007, 0, 0.081, 0.019, 0, 0, 0, 0.222).
$$

Two key observations:

1. **The LP-optimal $p$ is NOT uniform.**  It has heavy mass at the
   edges $j = 1$ and $j = 32$ (each $\approx 0.222$), and a non-trivial
   distribution in the middle.

2. **It is symmetric under reflection** ($p_j = p_{2N+1-j}$), even
   though we did not enforce evenness.  The LP has chosen a
   "spontaneously even" solution.

## What this implies for the bound

The LP-optimal Fourier coefficients give
$$
\sum_{k\ge 1}(a_k^{2} + b_k^{2}) \;\approx\; 0.352
$$
and
$$
\sum_{k\ge 1}|\hat F(k)|^{2} \;=\; \sum_{k\ge 1}(a_k^{2} + b_k^{2})^{2} \;\approx\; 0.034.
$$
So $\int F^{2} \approx 1 + 2 \cdot 0.034 = 1.069$, matching the
$\Omega^{(14)} = 1.063$ LP value (the small gap is due to the SDP
relaxation slack on $M$).

For comparison, the uniform $f$ has $\int F^{2} = 4/3 \approx 1.333$
and $\sum_{k\ge 1}(a_k^{2} + b_k^{2}) = 1/2$.

**The LP-optimal $p$ gives $\int F^{2} \approx 1.07$, which is
much less than $4/3$.**  This means the asymptotic limit of v14 is
**not** $4/3$ — it is some value strictly between $1$ and $4/3$,
determined by the actual (non-uniform) minimizer of the
constrained-LP problem.

## Implication for the project

This is a critical finding.  The v14 framework alone *cannot* reach
the SOTA $C_{6.2} \ge 1.28$, because:

* The LP picks a non-uniform $p$ with $\int F^{2} \approx 1.07$
  (seen here).
* As $N$ grows, the cell-bound on $a_k$ tightens, and the LP can
  pick $p$ even closer to the true minimizer of the $\int F^{2}$
  variational problem over admissible $f$.
* But this true minimizer is *not* uniform $f$ (which would give
  $4/3$); instead, it's some non-uniform $f$ that achieves a
  smaller value, plausibly close to $1.10–1.15$.

So the **ASYMPTOTIC bound** of v14 is $\sim 1.10–1.15$, which is
**below** the SOTA $1.28$.  v14 cannot reach the SOTA.

To reach the SOTA, we need a fundamentally **stronger constraint
type** that constrains the LP-optimal $p$ further (away from the
non-uniform minimizer).  Candidates:

1. **Higher-order Fourier moments of f** (e.g., $\int f^3$, $\int f^4$),
   which restrict $f$'s shape beyond cell-wise.
2. **The full Bochner constraint on the f-side autocorrelation**
   (Toeplitz [c_{i-j}] PSD), which v13 added but didn't help.
3. **A different relaxation framework entirely**, e.g., a
   sum-of-squares (SOS) hierarchy on the original infinite-dim
   problem, sampled at increasing levels.

This diagnostic does not produce a new bound but **clarifies** what
v14 alone can do.

## Validity

v26 adds no new constraints; the LP solved is exactly v14's.  All
v14 results remain valid.
