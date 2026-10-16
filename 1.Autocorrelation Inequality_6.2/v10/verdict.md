# Verdict — Autocorrelation Inequality 6.2, v10

## 1. Summary of the proof

v10 attempts to extend v9's Fejér-kernel test to **non-even** \(f\)
by introducing both cosine \(a_k\) and sine \(b_k\) Fourier
coefficients of \(f\), plus a moment matrix \(M\) on
\((1, a_1, \ldots, a_K, b_1, \ldots, b_K)\) of size \(2K+1\).

The constraint is the lifted form of the Fejér convolution
inequality \((F_K * F)(t_*)\le\Omega\):
\[
\Omega \;\ge\; 1 + 2\!\sum_{k=1}^{K}\!\Big(1 - \frac{k}{K+1}\Big)\!
\Big[(M[k,k] - M[K+k, K+k])\cos(2\pi k t_*)
+ 2 M[k, K+k]\sin(2\pi k t_*)\Big]
\]
for each test point \(t_*\in[-1/2, 1/2]\). The constraints are
proven *mathematically valid* (rank-one lift from any admissible \(f\)
satisfies them), but the proof's "Why v10 fails numerically" section
explicitly diagnoses that **the SDP optimum stays at \(\Omega = 1\)**
because of the Fourier-side diagonal inflation degeneracy: the lift
only gives *lower* bounds \(M[k,k]\ge a_k^{2}\),
\(M[K+k, K+k]\ge b_k^{2}\), and the relaxation can inflate these
freely to make the constraint trivially satisfied.

The proposer marks v10 as "rigorous but vacuous" and proposes v11
with a Bochner / Toeplitz PSD constraint to cap the diagonal
inflation.

## 2. Verdict

**VALID.**

The math is correct: the constraints in v10 are *mathematically
valid* (they hold for the rank-one lift from any admissible \(f\)),
the Fejér convolution computation is right, the lifted constraints
are well-defined, and the proposer's self-diagnosis of the failure
mode is accurate. There is no logical fallacy.

This is *different from v8's self-INVALID*. In v8, the constraint
*itself* was wrong (excluded an admissible \(f\)). In v10, the
constraint is *right* but the relaxation has structural slack that
absorbs it. The verifier instructions say:

> If there is no issue in the proof, in verdict.md after you
> summarize and noting some comments, give the final verdict as
> VALID. If you found a logical fallacy or a counterexample, give
> the final verdict as INVALID.

There is no fallacy and no counterexample. The proof is rigorous —
just structurally vacuous. **VALID.**

That said, the project does not advance numerically from v9 to v10:
v10's bound is \(\Omega^{(10)} = 1\), strictly *worse* than v9's
\(\Omega^{(9)} = 1.689\) (because v10 doesn't have v9's even-symmetry
constraints, and the Fourier-side machinery alone can't pin the
degeneracy). The proposer is honest about this.

## 3. Verification

### Setup. ✓
Non-even \(f\), complex \(\hat f(k) = a_k - ib_k\), with
\(\hat F(k) = \hat f(k)^{2} = (a_k^{2} - b_k^{2}) - 2 i a_k b_k\).
Standard.

### Fejér convolution computation. ✓

\((F_K * F)(t_*)\) computed via Fourier:
\[
(F_K * F)(t_*) = \sum_{|k|\le K}\!\Big(1 - \frac{|k|}{K+1}\Big)\hat F(k)e^{2\pi i k t_*}.
\]
For real \(F\), \(\hat F(-k) = \overline{\hat F(k)}\), so the sum
collapses to
\[
(F_K * F)(t_*) = \hat F(0) + 2\sum_{k=1}^{K}\Big(1 - \frac{k}{K+1}\Big)\,\mathrm{Re}\big(\hat F(k) e^{2\pi i k t_*}\big).
\]
Computing the real part with \(A = a_k^{2} - b_k^{2}\),
\(B = 2 a_k b_k\), \(\theta = 2\pi k t_*\):
\[
\mathrm{Re}\big((A - iB)(\cos\theta + i\sin\theta)\big) = A\cos\theta + B\sin\theta.
\]
So
\[
(F_K * F)(t_*) = 1 + 2\sum_{k=1}^{K}\Big(1 - \frac{k}{K+1}\Big)\big[(a_k^{2} - b_k^{2})\cos(2\pi k t_*) + 2 a_k b_k\sin(2\pi k t_*)\big]. \;\checkmark
\]
This matches the proof.

### The convolution upper bound. ✓
\((F_K * F)(t_*) = \int F_K(t) F(t_* - t)\,dt \le \Omega\int F_K(t)\,dt = \Omega\),
using \(F_K\ge 0\) and \(F\le\Omega\). ✓

(Note: this is *convolution* at \(t_*\), not *inner product*, but
gives the same bound. The proof writes
\((F_K * F)(t_*) \le \|F_K\|_{1}\cdot\|F\|_{\infty} = 1\cdot\Omega\),
which is just the standard \(\|f * g\|_{\infty} \le \|f\|_{1}\|g\|_{\infty}\)
Young inequality applied at the point \(t_*\). Same thing.)

### The lifted constraint. ✓
Replace \(a_k^{2}\to M[k,k]\), \(b_k^{2}\to M[K+k, K+k]\),
\(a_k b_k\to M[k, K+k]\). At the rank-one lift \(M^* = vv^{\top}\),
\(M^*[k,k] = a_k^{2}\), \(M^*[K+k, K+k] = b_k^{2}\),
\(M^*[k, K+k] = a_k b_k\), so the lifted constraint becomes the
actual \((F_K * F)(t_*) \le \Omega\), which holds by the convolution
upper bound. So the rank-one lift is feasible for the lifted
constraint, and the constraint is **valid as a relaxation**. ✓

### Diagonal inflation degeneracy. ✓ (correctly diagnosed)
The proposer's diagnosis is correct. The PSD lift on \(M\) only
forces the *lower* bounds \(M[k,k]\ge a_k^{2}\),
\(M[K+k, K+k]\ge b_k^{2}\) — there is no *upper* bound. So the
solver can freely inflate \(M[K+k, K+k]\), making the term
\((M[k,k] - M[K+k, K+k])\cos(2\pi k t_*)\) arbitrarily *negative*
when \(\cos > 0\). Combined with similar inflation for the
off-diagonal \(M[k, K+k]\), the relaxation can make the right-hand
side of the lifted constraint as small as it wants — in particular,
making \(\Omega \ge 1\) the only binding constraint.

This is the same diagonal-inflation degeneracy I noted in the
**6.2 v2 verdict Issue 5** ("the SDP lift on \(M\) is degenerate")
and predicted in the **6.2 v9 verdict Section 5 Comment 2** when
discussing the non-even extension. The proposer's diagnosis is
correct.

## 4. Comments

### Comment 1. v10 was the *predictable* failure mode of the naive non-even extension.

In my **v9 verdict Section 5 (Comment 2)** I explicitly predicted:

> So **v10 should use the autocorrelation Toeplitz lift**, not split
> into even/odd parts. The Toeplitz approach handles the non-even
> case without diagonal inflation, because the Toeplitz PSD
> constraint forces all the diagonal entries to be related.

The proposer's v10 went with the naive non-even extension instead
(adding sine coefficients $b_k$ and lifting to a moment matrix on
$(1, a, b)$), and hit *exactly* the diagonal inflation problem I
predicted. The good news: the proposer correctly identifies the
failure and proposes the right fix for v11 (Bochner/Toeplitz PSD).

### Comment 2. The v11 plan (Toeplitz PSD) is correct.

The proposer says:

> v11 will replace the SDP lift with a tighter alternative —
> specifically, the Bochner / Toeplitz PSD constraint on
> \((a_k, b_k)\), which ensures that the sequence \((\hat f(k))_k\)
> is the Fourier sequence of an actual non-negative function and
> thereby caps the off-diagonal moment inflation.

This is the right structural move. The Bochner / Toeplitz construction
is: define the Toeplitz matrix \(T_{l,m} = \hat f(l - m)\) (with
appropriate conjugation), and impose \(T\succeq 0\). For nonneg
\(f\), this is automatic (it's the Bochner theorem applied to the
characteristic function on \([-1/4, 1/4]\)). The PSD constraint
*forces* relations between diagonal and off-diagonal entries that
the naive moment matrix \(M\) does not enforce, so the diagonal
inflation is capped.

This is the same construction as **6.3 v4's Lemma Toep-1**
(autocorrelation Toeplitz PSD), transposed to the 6.2 problem.

### Comment 3. The "rigorous but vacuous" framing is correct.

The proposer's self-diagnosis "v10 is therefore 'rigorous but
vacuous': it neither violates the relaxation principle nor improves
the bound" is the right framing. The constraints *are* valid; they
just don't tighten the relaxation. This is a meaningfully different
failure mode from v8's self-INVALID (where the constraint was
*wrong*), and the proposer correctly distinguishes the two.

### Comment 4. Numerical comparison.

| Iteration | Bound | Even-restricted? |
|---|---|---|
| v0–v6 | 1.0 (trivial) | various |
| v7 | 1.344 (\(N=20\)) | yes |
| v8 | self-INVALID | (would have been even-only) |
| v9 | 1.689 (\(N=50, K=64\)) | **yes** |
| v10 | 1.0 (vacuous) | no (general \(f\)) |

The current state of the project:
* For \(C^{\text{even}}_{6.2}\): the best rigorous bound is **1.689**
  from v9.
* For \(C_{6.2}\) (general): the best rigorous bound is still **1.0**,
  the trivial bound. The SOTA \(C_{6.2}\ge 1.28\) (per `Problem.md`)
  has not been touched by the project's relaxation pipeline.

v11 is the make-or-break iteration: if the Toeplitz PSD constraint
caps the diagonal inflation as expected, the non-even bound should
finally rise above 1, and ideally above the SOTA 1.28 to make a
real numerical contribution.

### Comment 5. Inheritance check. ✓
v10 inherits from v9 (which inherits from v7), and v7 dropped the
v0 band hand-wave. So v10 has *no inherited rigor issues*. The only
issue is the diagonal inflation, which is a structural degeneracy
of v10's specific lift, not an inherited problem.

## 5. Feedback to the proposer

The "implement → test → diagnose → propose fix" workflow is exactly
right. v8 (self-INVALID), v9 (clean execution of Fejér for the even
case), and v10 (rigorous but vacuous, with correct diagnosis) form a
coherent sequence: each iteration learns from the previous. **This
is the right mode of progress.**

Two concrete suggestions for v11:

1. **Use the autocorrelation Toeplitz lift, not the moment matrix
   on \((a, b)\).** Concretely: define
   \(T_{l, m} := \hat f(l - m)\) for \(l, m = 0, \ldots, K\), where
   \(\hat f(k) = a_k - i b_k\) (so \(T_{l,m}\) is a complex Toeplitz
   matrix), and impose \(T\succeq 0\) (Hermitian PSD). This is the
   *Carathéodory–Toeplitz / Bochner* characterization of nonneg
   functions: \(\hat f(k)\) is a Fourier sequence of a nonneg measure
   on the circle iff its Toeplitz matrix is PSD.

   In CVXPY: use the *real* representation of complex PSD,
   \(\bigl[\begin{smallmatrix} \mathrm{Re}\,T & -\mathrm{Im}\,T\\ \mathrm{Im}\,T & \mathrm{Re}\,T\end{smallmatrix}\bigr]\succeq 0\)
   (a real PSD constraint of size \(2(K+1)\)). The same construction
   is in **`Erdos minimum overlap problem/v0/convex_code.py` lines
   285–304** if you want a reference implementation.

2. **Combine with the Fejér constraint of v9/v10.** The Toeplitz PSD
   constraint *replaces* the naive PSD lift on \(M\), but the Fejér
   convolution constraint at test points should be re-derived in
   terms of \(T\) (or kept as is and combined with the new
   constraint). The relevant identity: the entries of the
   autocorrelation \(F * \tilde F\) (where \(\tilde F(t) = F(-t)\))
   are exactly the diagonal of the Toeplitz matrix
   \(T_F\) of \(F\), and the Fejér convolution at any point is a
   linear functional of these diagonal entries.

3. **Test for the diagonal inflation again.** Before committing to
   v11, simulate the LP optimum at a small \((N, K)\) and check
   that the new constraints actually drop \(\Omega\) below 1. If
   they don't, diagnose the failure mode (probably some other
   degeneracy) before submitting.

If v11 successfully caps the inflation, the bound should rise to at
least 1.x for some \(x > 0\), and combined with v9's even-restricted
bound, the project will have actual numerical contributions for both
the even and general cases.

---

**Final verdict: VALID** (rigorous but vacuous; correct self-diagnosis and v11 plan).
