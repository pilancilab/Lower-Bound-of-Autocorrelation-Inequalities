# Verdict — Autocorrelation Inequality 6.2, v14

## 1. Summary of the proof

v14 finally executes the **autocorrelation Fejér construction** that
I recommended in the v9 verdict §5 (and re-recommended in v10, v12,
v13). The key idea: instead of testing the Fejér kernel against
\(F = f * f\) directly (which produces a constraint involving
\(a_k^{2} - b_k^{2}\), allowing diagonal inflation), test the Fejér
kernel against the **autocorrelation** \(F * \check F\) where
\(\check F(t) = F(-t)\). The autocorrelation has Fourier transform
\(|\hat F(k)|^{2} = (a_k^{2} + b_k^{2})^{2}\), with a *plus* sign that
makes the constraint immune to the diagonal-inflation degeneracy.

The constraint chain:
\[
\underbrace{\int F_K(F * \check F)\,dt}_{\text{Parseval} = 1 + 2\sum (1-k/(K+1))|\hat F(k)|^{2}}
\;\le\; \|F_K\|_{1}\cdot\|F*\check F\|_{\infty} = \int F^{2} \le \Omega\int F = \Omega.
\]
Lifting via two layers:
* \(u_k := M[k,k] + M[K+k, K+k] \ge a_k^{2} + b_k^{2}\).
* \(v_k \ge u_k^{2}\) (rotated SOC).

The lifted constraint is
\(\Omega \ge 1 + 2\sum_{k=1}^{K}(1 - k/(K+1))v_k\), linear in
\((\Omega, v_k)\) and convex via the SOC \(v_k \ge u_k^{2}\).

The proof concludes with numerical results showing
\(\Omega^{(14)}\) **strictly above 1** for every \((N, K)\) tested,
reaching **1.063 at \((N, K) = (16, 16)\)** — the *first*
non-trivial rigorous lower bound on the **unrestricted** \(C_{6.2}\)
in the project history.

## 2. Verdict

**VALID.**

This is the **third major milestone of the project** (after v7
dropping the v0 hand-wave for even \(f\), and v9 executing the
Fejér-kernel test for even \(f\)). v14 is the first iteration in
which:

1. The construction works for **all** admissible \(f\), not just
   even ones.
2. The numerical bound is **strictly above 1** without any of the
   inheritance issues from v0/v1/v2/v3/v4/v5 (v14 inherits only from
   v9 / v11, both of which are clean).
3. The diagonal-inflation degeneracy that defeated v10/v12/v13 is
   structurally bypassed by using \(|\hat F(k)|^{2}\) (which is a
   sum of squares, not a difference) instead of \(\hat F(k)\) (which
   has a signed real part).

The proof is rigorous, the math is correct, the lift is convex
(rotated SOC for \(v_k \ge u_k^{2}\)), and the numerical bound is
real. This is exactly the construction the v9 verdict §5 pointed
to:

> **The right fix** is to use a *positive* linear combination of
> coefficient terms — i.e., to use the *autocorrelation*
> \(F * \tilde F\) where \(\tilde F(t) = F(-t)\), which is
> automatically nonneg in Fourier space (\(|\hat F(k)|^{2}\)) and
> gives a Fejér-style constraint of the form
> \(\Omega^{2} \ge \sum_{|k|\le K}(1 - |k|/(K+1))|\hat F(k)|^{2}\).

(I had it as \(\Omega^{2}\) on the LHS, expecting the bound to
involve squaring, but v14's actual derivation gets just \(\Omega\)
on the LHS via the cleaner chain
\(\int F^{2} \le \Omega\int F = \Omega\). v14's version is tighter
than my prediction.)

## 3. Verification

### Property (A): \(\widehat{F * \check F}(k) = |\hat F(k)|^{2}\). ✓
Convolution theorem: \(\widehat{F*\check F}(k) = \hat F(k)\hat{\check F}(k)\).
Computing \(\hat{\check F}(k) = \int F(-t)e^{-2\pi i k t}\,dt = \int F(s)e^{2\pi i k s}\,ds = \overline{\hat F(k)}\)
(via \(s = -t\)). So
\(\widehat{F*\check F}(k) = \hat F(k)\overline{\hat F(k)} = |\hat F(k)|^{2}\),
real and nonneg. ✓

### Property (B): \(F * \check F \ge 0\). ✓
Both \(F \ge 0\) and \(\check F \ge 0\), so their convolution is nonneg. ✓

### Property (C): \((F * \check F)(0) = \int F^{2}\) and the autocorrelation peaks at 0. ✓

\((F * \check F)(0) = \int F(s)\check F(-s)\,ds = \int F(s) F(s)\,ds = \int F^{2}\). ✓

The peak claim: by Cauchy–Schwarz,
\(|(F * \check F)(t)| = |\int F(s) F(s+t)\,ds| \le \big(\int F(s)^{2}\,ds\big)^{1/2}\big(\int F(s+t)^{2}\,ds\big)^{1/2} = \int F^{2} = (F*\check F)(0)\),
where the second equality uses translation invariance of \(L^{2}\) norm. ✓

### The constraint chain. ✓
\[
\int F_K\,(F * \check F)\,dt \;\le\; \|F_K\|_{1}\cdot\|F*\check F\|_{\infty} = 1\cdot\int F^{2} = \int F^{2} \;\le\; \Omega\int F = \Omega\cdot 1 = \Omega.
\]
* First inequality: Hölder \(\int gh \le \|g\|_{1}\|h\|_{\infty}\) with both \(g = F_K, h = F*\check F\) nonneg, and \(\|F_K\|_{1} = 1\) (Property F2 from v9).
* Second equality: Property (C).
* Third inequality: \(F \ge 0\), \(F \le \Omega\), so \(\int F\cdot F \le \int \Omega \cdot F = \Omega\int F = \Omega\).

All correct. ✓

### Parseval expansion of the LHS. ✓
\[
\int F_K(F*\check F)\,dt = \sum_{|k|\le K}\widehat{F_K}(k)\overline{\widehat{(F*\check F)}(k)} = \sum_{|k|\le K}\Big(1 - \frac{|k|}{K+1}\Big)|\hat F(k)|^{2},
\]
both \(\widehat{F_K}\) and \(\widehat{F*\check F}\) being real. The
sum collapses by symmetry (\(|\hat F(-k)|^{2} = |\hat F(k)|^{2}\) and
\(|\hat F(0)|^{2} = (\int F)^{2} = 1\)) to:
\[
1 + 2\sum_{k=1}^{K}\Big(1 - \frac{k}{K+1}\Big)|\hat F(k)|^{2}. \;\checkmark
\]

### The constraint (★): \(1 + 2\sum_{k=1}^{K}(1 - k/(K+1))|\hat F(k)|^{2} \le \Omega\). ✓
By chaining the previous two computations.

### The lift (★★): all variables nonneg, no minus sign. ✓

The crucial computation: \(|\hat F(k)|^{2} = |\hat f(k)|^{4} = (a_k^{2} + b_k^{2})^{2}\),
where the **plus** sign comes from the *modulus* (not the real
part). At every step in the lift, the variables are nonneg sums:
\(u_k = M[k,k] + M[K+k, K+k] \ge a_k^{2} + b_k^{2}\), and
\(v_k \ge u_k^{2}\). 

**This is the key structural difference from v10/v12/v13:** in those
versions, the LP saw \(\hat F(k)\) directly (which has both real and
imaginary parts), and the real part \(\mathrm{Re}\,\hat F(k) = a_k^{2} - b_k^{2}\)
has a minus sign that the LP could exploit by inflating
\(M[K+k, K+k]\). In v14, the LP sees \(|\hat F(k)|^{2} = (a_k^{2} + b_k^{2})^{2}\),
which is *intrinsically nonneg* and *monotone increasing* in
\(M[K+k, K+k]\), so inflation only *hurts* the LP (makes \(\Omega\)
larger), not helps.

### Convexity of the lift. ✓
* \(u_k = M[k,k] + M[K+k, K+k]\): linear in \(M\).
* \(v_k \ge u_k^{2}\): rotated SOC; CVXPY DCP atom `cp.square` accepts this.
* Constraint \(\Omega \ge 1 + 2\sum (1 - k/(K+1))v_k\): linear in \((\Omega, v_k)\).
* PSD lift on \(M\): SDP constraint.
* f-side cell-bound on \((a_k, b_k)\) from \(p_j\): linear (inherited from v9/v12).

All constraints are convex. The full program is a convex SDP. ✓

### Validity of the rank-one lift. ✓
For any admissible \(f\) (not just even):
* \(p_j^{*} = \int_{I_j} f\), \(a_k^{*} = \int f\cos\), \(b_k^{*} = \int f\sin\).
* \(M^{*} = v v^{\top}\) with \(v = (1, a_1^{*}, \ldots, a_K^{*}, b_1^{*}, \ldots, b_K^{*})\).
* \(u_k^{*} = (a_k^{*})^{2} + (b_k^{*})^{2} = |\hat f(k)|^{2}\), at the PSD-lift boundary.
* \(v_k^{*} = (u_k^{*})^{2} = |\hat f(k)|^{4} = |\hat F(k)|^{2}\), at the SOC boundary.
* Constraint (★★) reduces to (★) at the rank-one lift, which holds
  by the derivation of (★) above. So
  \(\Omega^{*} = \|f * f\|_\infty\) satisfies the constraint, and the
  rank-one lift is feasible.

Therefore \(\Omega^{(14)} \le \|f * f\|_\infty\) for every admissible
\(f\), giving \(\Omega^{(14)} \le C_{6.2}\). ✓

## 4. Numerical results

| \(N\) | \(K\) | \(\Omega^{(14)}\) |
|---|---|---|
| 8 | 4 | 1.026 |
| 8 | 8 | 1.033 |
| 8 | 16 | 1.036 |
| 16 | 4 | 1.040 |
| 16 | 8 | 1.055 |
| 16 | 16 | **1.063** |

**Strictly above 1 for every \(N \ge 8\), \(K \ge 4\) tested.** This
is the first non-even rigorous bound on \(C_{6.2}\) in the project
history. Comparison with prior iterations:

| Iteration | Best bound | Even-restricted? |
|---|---|---|
| v0–v6 | 1.0 (trivial) | various |
| v7 | 1.344 (\(N=20\)) | yes |
| v9 | 1.689 (\(N=50, K=64\)) | yes |
| v11 | 1.376 (\(N=8, K=16\)) | yes |
| v10/v12/v13 | 1.0 (vacuous) | no |
| **v14** | **1.063 (\(N=16, K=16\))** | **no** |

The v14 bound 1.063 is well below the SOTA \(C_{6.2}\ge 1.28\) at
the small \((N, K)\) tested. The proposer notes that larger
\((N, K)\) runs are needed to see the asymptotic limit, and v15
will explore larger parameters and add the v11 windowed (W)
constraints to tighten further.

## 5. Comments

### Comment 1. The autocorrelation lift is the right structural fix.

I've been recommending this since the **v9 verdict §5 Comment 2**:

> **The right fix** is to use a *positive* linear combination of
> coefficient terms — i.e., to use the *autocorrelation* \(F * \tilde F\)
> where \(\tilde F(t) = F(-t)\), which is automatically nonneg in
> Fourier space (\(|\hat F(k)|^{2}\)).

v10, v12, v13 all attempted variations on the naive non-even Fourier
extension and hit the diagonal inflation. v14 finally executes the
autocorrelation construction, and the diagonal inflation is
structurally bypassed.

The key insight: \(|\hat F(k)|^{2} = (a_k^{2} + b_k^{2})^{2}\) is
*intrinsically nonneg and monotone increasing in both \(a_k^{2}\) and
\(b_k^{2}\)*. Inflating \(M[K+k, K+k]\) (which represents \(b_k^{2}\)
in the lift) only makes the constraint *tighter*, never looser. So
the LP cannot escape the constraint by inflating; the optimum sits
at the rank-one boundary.

### Comment 2. The numerical bound is small but real.

\(\Omega^{(14)} = 1.063\) at \((N, K) = (16, 16)\). The proposer
correctly notes that this is well below the SOTA \(C_{6.2}\ge 1.28\)
and that larger \((N, K)\) runs are needed.

A few things that should help v14's bound grow at scale:

1. **Larger \(K\)**: the Fejér constraint at \(K = 16\) only sees
   the *first 16* Fourier coefficients of \(F\). Going to \(K = 32, 64, 128\)
   should significantly tighten the bound, since the autocorrelation
   uses *all* the coefficients.
2. **Larger \(N\)**: the f-side cell discretization at \(N = 16\) is
   coarse. v9's best (1.689) was at \(N = 50\), so v14 should
   easily go to \(N = 30, 50, 100\).
3. **Combination with the v11 windowed (W) constraints**: v11
   tightened v9 by 0.027 with shifted windows; combining with v14
   should give similar tightening.
4. **The v13 Bochner constraint** also fits cleanly into v14's
   framework (it's a constraint on the first row of \(M\), and
   v14's lift uses the diagonal of \(M\) plus the SOC; they don't
   conflict). Adding both should help.

I expect v15/v16 with all of these combined to reach a bound in the
range 1.2–1.5 — possibly *meeting or exceeding the SOTA* 1.28 if
the convergence is favourable.

### Comment 3. The convexity is preserved by the rotated SOC.

The constraint \(v_k \ge u_k^{2}\) is a *rotated second-order cone*
constraint (equivalent to
\(\bigl[\begin{smallmatrix} v_k & u_k \\ u_k & 1\end{smallmatrix}\bigr]\succeq 0\)
or to the cone \(\{(v, u, 1) : v\cdot 1 \ge u^{2}\}\)). CVXPY's
`cp.square` atom accepts this directly. The full program is therefore
a convex SDP+SOC, solvable by MOSEK or similar.

### Comment 4. v14 vs v9: different metrics, same machinery.

v9 (the previous milestone) used the Fejér kernel test on \(F\)
directly, restricted to even \(f\). It got \(\Omega^{(9)} = 1.689\)
at \((50, 64)\) on the *even-restricted* constant \(C_{6.2}^{\text{even}}\).

v14 uses the Fejér kernel test on the *autocorrelation* \(F * \check F\),
with no even restriction. It gets \(\Omega^{(14)} = 1.063\) at
\((16, 16)\) on the *unrestricted* \(C_{6.2}\).

The two are not directly comparable: v9 is bounding a *larger*
constant (since \(C_{6.2}^{\text{even}} \ge C_{6.2}\)) and gets a
larger bound; v14 is bounding the *smaller* unrestricted constant
and gets a smaller bound. Both are correct — they just measure
different things.

For the project's main numerical target (the unrestricted
\(C_{6.2}\)), v14 is the relevant bound. v9 is informational about
the even-restricted case but doesn't bound the SOTA target.

### Comment 5. The autocorrelation derivation is *cleaner* than my v9 prediction.

In the v9 verdict §5 Comment 2, I predicted the autocorrelation
constraint would have the form:
\[
\Omega^{2}\ge \sum_{|k|\le K}(1 - |k|/(K+1))|\hat F(k)|^{2}
\]
(with \(\Omega^{2}\) on the LHS, expecting the bound to come from
\(\int(F*\check F) = (\int F)^{2}\cdot |\text{stuff}|\) of degree 2).

v14's actual derivation is *cleaner*: the LHS is just \(\Omega\)
(not \(\Omega^{2}\)), via the chain
\(\int F^{2} \le \Omega\int F = \Omega\). So v14's constraint is
*linear in \(\Omega\)* on the LHS, which fits cleanly into the
LP/SDP framework without needing a Schur complement on \(\Omega\)
itself.

The proposer's improvement on my prediction is the right move.

## 6. Feedback to the proposer

This is the most significant iteration in 6.2 since v9. The
autocorrelation Fejér construction is the right structural fix for
the non-even case, and v14 executes it cleanly. The math is rigorous,
the numerical bound is non-trivial, and the inheritance is clean.

Three concrete suggestions for v15:

1. **Larger parameters.** Run v14 at \(N \in \{32, 64, 100\}\) and
   \(K \in \{32, 64, 128\}\). The bound should grow significantly,
   possibly approaching or exceeding the SOTA 1.28. This is the most
   important next step — *purely a parameter sweep*, no new
   theoretical work.
2. **Combine with v11's shifted windows.** The (W) family from v11
   is independent of the autocorrelation Fejér constraint and adds
   strict information. Adding both should give a larger bound.
3. **Combine with v13's Bochner Toeplitz.** Same: the Bochner
   constraint on \(\hat f\) is independent of the autocorrelation
   constraint on \(|\hat F|^{2}\); adding both should help.

If v15 does (1) + (2) + (3) and the bound exceeds 1.28, the project
will have a real numerical contribution to the SOTA on \(C_{6.2}\).

A small note on terminology: the autocorrelation construction here
is essentially the "Wiener–Khinchin" identity applied to \(F\). It's
worth noting in the proof or in v15 that the same structure underlies
the **6.3 v4 Toeplitz Lemma** (where the autocorrelation Toeplitz
matrix is forced PSD by Bochner's theorem). The two problems share
the same structural ingredient; v14 brings this insight to 6.2.

---

**Final verdict: VALID.** This is the **third major milestone** of the project.
