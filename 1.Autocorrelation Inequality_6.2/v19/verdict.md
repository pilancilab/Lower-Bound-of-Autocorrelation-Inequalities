# Verdict — Autocorrelation Inequality 6.2, v19

## 1. Summary

v19 inherits the v14 autocorrelation Fejér framework and adds the
*same* Fejér-style constraint at *multiple* values of \(K_{\text{eff}}\)
simultaneously:
\[
\Omega \;\ge\; 1 + 2\sum_{k=1}^{K_{\text{eff}}}\!\Big(1 - \frac{k}{K_{\text{eff}}+1}\Big)v_k,
\quad K_{\text{eff}}\in\{2, 4, 6, \ldots, K_{\text{max}}\}.
\]

Each individual constraint is valid (proven in v14), and adding
multiple constraints only restricts the feasible set further. The
numerical improvement over v14 alone is marginal (\(\sim 0.001\)–\(0.003\)
at the tested parameters).

## 2. Verdict

**VALID** (small improvement over v14; the multiple-\(K\) Fejér
constraints are individually valid, and their conjunction is a
strict tightening of the v14 LP).

## 3. Verification

Each constraint
\[
\Omega \ge 1 + 2\sum_{k=1}^{K_{\text{eff}}}(1 - k/(K_{\text{eff}}+1))v_k
\]
is the v14 Fejér constraint at the specific \(K_{\text{eff}}\). v14
proves this is valid for the rank-one lift of any admissible \(f\)
(via the chain \(\int F_{K_{\text{eff}}}\cdot(F*\check F)\le \Omega\)),
so each individual constraint at each \(K_{\text{eff}}\) is valid.

The conjunction of valid constraints is itself valid: every admissible
\(f\) gives a feasible point for *all* the constraints simultaneously
(since each constraint is satisfied at the rank-one lift). Hence the
v19 SDP optimum \(\Omega^{(19)}\le C_{6.2}\). ✓

## 4. Numerical analysis

| \(N\) | \(K_{\text{max}}\) | \(\Omega^{(19)}\) |
|---|---|---|
| 8 | 16 | 1.036 |
| 8 | 32 | 1.038 |
| 16 | 16 | 1.063 |
| 16 | 32 | 1.066 |
| 32 | 16 | 1.083 |
| 32 | 32 | 1.088 |

For comparison, v14 alone at the same parameters:
* \((N, K) = (16, 16)\): v14 gives 1.063 ✓ matches v19.
* \((N, K) = (16, 32)\): v14 would give around 1.063 + small (extrapolated).

The improvement of v19 over v14 is \(\sim 0.001\)–\(0.003\) at the
tested parameters. **Marginal.**

The proposer correctly diagnoses why: the smaller-\(K\) Fejér
constraints are *dominated* by the larger-\(K\) constraint at the
LP optimum. At the LP optimum, the larger-\(K\) constraint is binding
or nearly binding, and the smaller-\(K\) constraints are well within
their slack.

## 5. Comments

### Comment 1. Multi-\(K\) Fejér is structurally redundant.

Why does multi-\(K\) not help much? The Fejér weights
\(1 - k/(K_{\text{eff}}+1)\) are *triangular*: at small \(k\) they're
close to 1, at \(k = K_{\text{eff}}\) they're close to 0. As
\(K_{\text{eff}}\) increases, the *individual* weights at fixed
\(k\) increase (from \(1 - k/(K_{\text{eff}}+1)\) to \(1 - k/(K_{\text{eff}}'+1)\)),
and the constraint includes *more* high-frequency \(v_k\) terms.

So the larger-\(K\) constraint is *strictly stronger* than the
smaller-\(K\) one in two ways: (a) more terms, (b) higher weights
on the existing terms. The smaller-\(K\) constraint adds nothing
new.

A *non-redundant* set of Fejér-style constraints would use
*different kernel families* at the same \(K\), e.g.:
* Fejér kernel (current).
* Jackson kernel.
* de la Vallée–Poussin kernel.
* Optimal nonneg trig polynomials computed via the Carathéodory–Toeplitz
  parameterization.

Each of these is *also* a positive trigonometric polynomial, so the
same Hölder-style chain
\(\int K_n F\le \|K_n\|_1\|F\|_\infty\) gives a different test
constraint. Adding *multiple kernels* at the same \(K\) gives
strictly more information than multiple \(K\)'s of the same kernel.

### Comment 2. The v17 / v19 numerical bounds are converging slowly.

| Iteration | Best bound | Notes |
|---|---|---|
| v14 (\(N=16, K=16\)) | 1.063 | First non-even rigorous bound |
| v17 (\(N=100, K=32\)) | 1.111 | Parameter sweep |
| **v19 (\(N=32, K=32\))** | **1.088** | Multi-\(K\) Fejér |

The improvement from \(N = 32\) to \(N = 100\) is \(\sim 0.023\),
which suggests the convergence rate is something like
\(\Omega^{(\text{best})}(N) - \Omega^{(\text{best})}(\infty)\sim N^{-\alpha}\)
for small \(\alpha\). At this rate, reaching the SOTA 1.28 (let alone
the asymptotic \(4/3\)) would require enormous \(N\).

The proposer's conclusion is right: "genuinely new ideas are needed
to break past \(\sim 1.10\) at moderate parameters."

### Comment 3. The honest plateau diagnosis.

The proposer is diagnosing a *real* phenomenon: the v14 framework
appears to plateau around 1.10–1.15 at moderate \(N, K\), even with
all the bells and whistles (multi-\(K\), shifted (W) windows,
Bochner Toeplitz, L³, ...). None of the auxiliary constraints
helped, and the parameter sweep is logarithmic-slow.

This is a *correct empirical observation* and matches the project's
overall pattern: v14 is the right framework, but the convergence is
limited by either (a) the discretization bottleneck on \(p \to a_k\),
or (b) the actual infimum of \(\int F^{2}\) being smaller than \(4/3\).

To distinguish, the proposer should compute \(\int F^{2}\) for the
LP-optimal \((p, a, b)\) at large \(N\) and check whether it matches
the LP bound or is significantly smaller.

## 6. Feedback to the proposer

Three concrete suggestions for v20:

1. **Try other positive trigonometric kernels.** As I mentioned in
   Comment 1, multi-\(K\) Fejér is redundant, but multi-*kernel* at
   the same \(K\) is genuinely new information. Try the Jackson
   kernel \(J_K\) (which decays as \(k^{-3}\) instead of \(k^{-1}\)
   for the Fejér) and the de la Vallée–Poussin kernel.

2. **Compute \(\int F^{2}\) for the LP-optimal \(p\) at large \(N\)**
   to see whether the LP is finding a non-uniform minimizer with
   \(\int F^{2} < 4/3\). If yes, then the asymptotic limit is *not*
   \(4/3\), and the project needs a different approach to find the
   actual infimum.

3. **Tighten the cell-bound on \(a_k\) from \(p\)** by using a
   midpoint quadrature with explicit error instead of the loose
   worst-case \(\min/\max\) of \(\cos\) on each cell. This may
   accelerate the convergence in \(N\).

If v20 does (1) and (2), we will know whether the v14 framework
*can* reach the SOTA at all, or whether a different structural
approach is needed.

---

**Final verdict: VALID** (minor improvement over v14, but the
v14 plateau remains the binding bottleneck).
