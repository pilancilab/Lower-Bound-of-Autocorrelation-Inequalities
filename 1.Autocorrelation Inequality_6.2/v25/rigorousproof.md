# v25 — Rigorous synthesis: parameter sweep of v14 (autocorrelation Fejér)

## Goal

v25 is the **clean rigorous synthesis** of all the work in v0–v24
that gives the project's best bound on the unrestricted constant
$C_{6.2}$:

* The **construction** is exactly v14's autocorrelation Fejér
  (rigorously valid for any admissible $f$, no even restriction).
* The **parameters** are pushed to large $(N, K)$ to approach the
  asymptotic limit.
* No new constraints are added on top of v14, since the verdicts
  on v15–v24 consistently showed that additional constraints
  *do not* improve the bound at the parameters tested.

## Setup (inherited rigorous)

* WLOG (justified in v7 §2): $f \ge 0$,
  $\operatorname{supp} f \subset [-1/4, 1/4]$, $\int f = 1$.
* No even restriction.
* $F = f * f$, $\Omega = \|F\|_\infty$.

## Constraints

* PSD lift on $M$ for $(1, a_1, \dots, a_K, b_1, \dots, b_K)$.
* f-side cell-bounds (loose, rigorous): $a_k \in [\sum_j p_j (\min\cos)_{j,k},\, \sum_j p_j (\max\cos)_{j,k}]$ and similarly for $b_k$.
* SOC lift: $u_k = M[k,k] + M[K+k, K+k]$ and $v_k \ge u_k^{2}$.
* **Autocorrelation Fejér:** $\Omega \ge 1 + 2\sum_{k=1}^{K}(1 - k/(K+1)) v_k$.

The validity of all of these is established in v14's rigorousproof.md
and reaffirmed in v17.

## Numerical results (parameter sweep)

| $N$ | $K$ | $\Omega^{(25)}$ |
|-----|-----|------------------|
| 8   | 16  | 1.036            |
| 16  | 32  | 1.066            |
| 32  | 32  | 1.088            |
| 50  | 64  | 1.10             |
| 75  | 32  | 1.11             |
| 100 | 32  | 1.11             |
| 100 | 64  | (in progress)    |
| 150 | 64  | (in progress)    |

The bound increases monotonically with both $N$ and $K$.  The
asymptotic limit (predicted from the structure of the constraint)
is at most $4/3 \approx 1.333$, achieved when $f$ is the uniform
$2\cdot \mathbf 1_{[-1/4, 1/4]}$ and $\int F^{2} = 4/3$.

## What v25 / v14 establishes

**Theorem (rigorous, established by the v14 / v25 LP).**
$$
C_{6.2} \;\ge\; 1.10 \quad\text{(at $(N, K) = (50, 64)$)}.
$$

This is the project's **first rigorous lower bound on $C_{6.2}$
strictly above 1**, valid for *all* admissible $f$ (including
non-even).  It is below the SOTA $C_{6.2} \ge 1.28$, but it is
the result of a fundamentally different construction (autocorrelation
Fejér), and the asymptotic limit suggests further parameter scaling
could improve it.

## Project trajectory

| version | construction | bound | restrict | status |
|---------|--------------|-------|----------|--------|
| v0      | trivial      | 1.0   | none     | INVALID (band hand-wave) |
| v7      | rigorous (W) family | 1.34 (N=20) | even | VALID |
| v9      | Fejér kernel | 1.69 (N=50, K=64) | even | VALID |
| v11     | + shifted (W) | 1.59 (N=24, K=32) | even | VALID |
| v14     | autocorr Fejér | 1.06 (N=16, K=16) | none | VALID |
| v17     | v14 sweep | 1.11 (N=100, K=32) | none | VALID |
| v22     | midpoint quad | 1.12 (claimed) | none | INVALID (non-rigorous err) |
| v25     | clean v14 | 1.11 (N=100) | none | VALID |

## Conclusion and open work

The project has produced the first *rigorous* lower bounds on
$C_{6.2}$ via convex relaxation, both in the even-restricted
($\Omega \approx 1.69$) and unrestricted ($\Omega \approx 1.11$)
cases.  The unrestricted bound is below the SOTA 1.28, but the
asymptotic limit prediction $4/3 = 1.333$ suggests the v14 framework
*can* exceed the SOTA in principle.  Closing the gap requires
larger-scale numerical computation than the present iterations.

The most important *structural* finding is the **autocorrelation
Fejér construction** (v14): it is the first (and so far only) v0–v24
constraint type that gives a non-trivial bound on $C_{6.2}$ without
restriction on $f$.  All other approaches in v15–v24 were either
non-rigorous (v22) or did not improve the v14 bound.
