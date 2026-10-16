# Verdict — Autocorrelation Inequality 6.2, v20

## 1. Summary

v20 adds a **Lasserre level-2 lift** on the \((a, b)\) variables: a
moment matrix \(M_2\) of size \(\binom{2K+1}{2}+2K+1\) indexed by
multi-indices of degree \(\le 2\), constrained \(M_2\succeq 0\).
Standard Lasserre level-2 SDP machinery.

The numerical bound is *identical* to v14 at the same parameters
(1.026 at \((8, 4)\), 1.040 at \((16, 4)\)). The level-2 lift does
not improve.

## 2. Verdict

**VALID** (rigorous level-2 Lasserre lift, but no numerical
improvement over v14).

The construction is standard:
* The rank-one moment sequence \(y_{\alpha} = a^{\alpha_a}b^{\alpha_b}\)
  for \(|\alpha|\le 2\) gives a rank-one PSD moment matrix
  \(M_2\succeq 0\) at the rank-one truth.
* In the relaxation, drop the rank-one requirement and impose
  \(M_2\succeq 0\) only.

This is a *valid* relaxation tightening: every admissible \(f\)
gives a feasible point of v20, so \(\Omega^{(20)}\le C_{6.2}\).

## 3. Why v20 doesn't help

The proposer correctly identifies the issue:
> The LP optimum at v14 already corresponds to a rank-one moment
> matrix on \((1, a, b)\), and lifting to higher-degree moments
> simply adds more rank-one PSD constraints that are automatically
> satisfied.

This is the standard "Lasserre tightness at low degree" phenomenon.
If the LP optimum is already at a rank-one moment matrix at level 1,
then the level-2 (and higher) Lasserre constraints are *automatically
satisfied* at that point, and the relaxation doesn't tighten.

To get level-2 to bind, you'd need a constraint that the rank-one
level-1 moment matrix *violates* but the rank-one level-2 *also*
satisfies — which is impossible if the v14 LP is optimal at a
rank-one truth.

## 4. Comments

### Comment 1. v20 is in the same equivalence class as v15, v16, v18, v19, v20.

The pattern is now overwhelming:

| Iteration | What it adds | Improvement over v14 |
|---|---|---|
| v15 | v7/v11 (W) windows | 0 |
| v16 | v13 Bochner Toeplitz | 0 |
| v18 | L³ constraint | 0 |
| v19 | Multi-\(K\) Fejér | ~0.003 |
| **v20** | **Lasserre level-2 on \((a, b)\)** | **0** |

**Each additional constraint is valid but structurally vacuous on
top of v14.** The v14 autocorrelation Fejér is the binding
constraint, and nothing else added so far moves the needle.

This is informative: it tells the proposer that the bottleneck is
*not* the discretization or the SDP lift, but rather the *fundamental
structural strength* of the v14 constraint at the working
parameters. To break past 1.10, the proposer needs a *qualitatively
different* constraint, not a refinement of the existing ones.

### Comment 2. What might actually help.

I count the project's "binding bottleneck" as the *infimum*
\(\inf_f \int(f*f)^{2}\) where the infimum is over admissible \(f\).
The proposer's prediction (uniform \(f\) achieves \(4/3\)) is
sound, but the LP convergence to this asymptotic limit is slow.

The reason for slow convergence may be:
1. **Discretization bottleneck**: the cell-bound on \(a_k\) from
   \(p\) is loose for large \(k\), so finer grids help only slowly.
2. **Non-uniform optimum**: the LP finds a non-uniform \(p\) with
   \(\int F^{2} < 4/3\), so the asymptotic limit is *not* \(4/3\)
   but something smaller.

To distinguish: extract the LP-optimal \((p, a, b)\) at large \(N\)
and compute \(\int F^{2}\) numerically (or compute it from the SDP
dual). If it matches the LP value, then (2) is the issue. If not,
(1) is the issue.

Either way, the path forward is *not* "add more SDP constraints".
It's either:
* Tighten the cell-bound on \(a_k\) via midpoint quadrature with
  explicit error.
* Or analytically determine \(\inf_f \int(f*f)^{2}\) (which is a
  separate variational problem) and check whether the LP is
  approaching it.

### Comment 3. Suggestion for v21.
**Stop adding constraints.** The next iteration should focus on:
1. Diagnosing the convergence bottleneck (which of the two issues
   above).
2. Either tightening the cell-bound (option 1) or analytically
   computing the infimum (option 2).

If the proposer continues adding constraints in the v15-v20 mode,
they'll just keep getting "VALID, no improvement" for each one.

---

**Final verdict: VALID** (no improvement over v14).
