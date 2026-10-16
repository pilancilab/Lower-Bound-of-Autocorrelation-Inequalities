# Verdict — Autocorrelation Inequality 6.2, v2

## 1. Summary of the proof

v2 inherits the v1 SDP (PSD lift on \(Q = pp^\top\), marginals,
nonnegativity, band constraints) and adds Fourier-side machinery.
Specifically:

* **Setup.** Period-1 Fourier coefficients of \(f\) on
  \([-1/4,1/4]\subset[-1/2,1/2]\): \(a_k,b_k\) with \(a_0=1,b_0=0\),
  and the identity \(\hat F(k) = \hat f(k)^2\) where \(F = f*f\). The
  expansion of \(F(t)\) on \([-1/2,1/2]\) in real cosines/sines (★) is
  derived from this.

* **(D) Cell-wise interval bounds on \(a_k, b_k\).** With
  \(\alpha^\pm_{j,k}=\max/\min_{x\in I_j}\cos(2\pi k x)\) and
  \(\beta^\pm_{j,k}\) similarly for \(\sin\):
  \[
  \sum_j p_j\,\alpha^-_{j,k}\le a_k\le\sum_j p_j\,\alpha^+_{j,k},
  \qquad
  \sum_j p_j\,\beta^-_{j,k}\le b_k\le\sum_j p_j\,\beta^+_{j,k}.
  \]
  This is the same template as White's Lemma 5 transposed to the
  problem 6.2 setting.

* **(E) Moment-matrix PSD lift.** \(M = vv^\top\) for
  \(v=(1,a_1,\dots,a_K,b_1,\dots,b_K)\); drop rank-one in the
  relaxation.

* **(F) Test-point upper bounds via the lifted \(F\).** For test points
  \(t_*\in[-1/2,1/2]\), the proof writes a partial-sum inequality
  asserting that the \(K\)-th symmetric partial sum of \(F\) at
  \(t_*\), evaluated through \(M\), is \(\le \Omega\).

## 2. Verdict

**INVALID.**

(D) and (E) are correct and rigorous. (F) has two independent
problems, one of which the proof itself acknowledges:

* **(F)'s naive form is false** because Fourier partial sums of \(F\)
  can overshoot \(F\) (Gibbs phenomenon), so partial-sum \(\le F(t_*)\)
  is **not** a universal inequality.
* The proof's *repaired* form requires a tail bound \(T_K\) that in
  turn needs an a-priori upper bound \(M_2\) on \(\int f^2\), and the
  LP supplies no such bound. The proof leaves \(M_2\) as a "known"
  external parameter without saying what it is.
* The proof then explicitly disclaims the rigorous form: "in the
  present implementation we ... simply add the lifted partial-sum
  inequality without the explicit tail slack. **Under that
  simplification the relaxation may slightly violate mathematical
  rigor by an amount bounded by the truncation tail.**" This is
  exactly the kind of "approximation that may or may not be a
  mathematically valid lower bound" the verifier instructions
  explicitly forbid.

In addition, v2 inherits without addressing the unresolved Issue 1
from the v1 verdict (the false band-integral identity in the
inherited band constraints).

## 3. Issues

### Issue 1 (FATAL). (F) as actually implemented is admittedly non-rigorous.

The proof writes (F) as
\[
\Omega \;\ge\; 1 + 2\sum_{k=1}^{K}\big[(M[k,k]-M[K+k,K+k])\cos(2\pi k t_*) + 2M[k,K+k]\sin(2\pi k t_*)\big]
\]
and identifies the right-hand side as the \(K\)-th symmetric partial
sum of the Fourier series of \(F\) at \(t_*\). It then asserts this
partial sum is \(\le F(t_*)\le \Omega\).

**The first inequality is false in general.** Fourier partial sums of
piecewise-smooth functions overshoot at jumps and corners (Gibbs
phenomenon). For nonneg \(F = f*f\) of any non-symmetric \(f\), the
partial sum can exceed \(F(t_*)\) at the kinks. Concretely: take any
of the \(f\)'s in the v1 verdict counterexample
(\(f = 4\cdot \mathbf 1_{[-1/8, 1/8]}\) on the \(N=1\) grid; \(F\) is a
triangle on \([-1/4,1/4]\) with peak height \(4\) at \(0\) and corners
at \(\pm 1/4, 0\)). The Fourier partial sums of this triangle at the
*corners* exceed the triangle's value by up to the Gibbs constant
\(\sim 0.09\) of the jump in the derivative. So
"partial sum \(\le F(t_*)\)" fails at \(t_*\in\{-1/4, 0, 1/4\}\).

If you then use these \(t_*\)'s as test points in (F), the constraint
\(\Omega \ge \) partial sum is *stronger* than \(\Omega\ge F(t_*)\),
i.e. it can rule out the rank-one lift \(M=vv^\top\) coming from the
true \(f\) — meaning the LP feasible set fails to contain the
discretized image of every admissible \(f\), and the LP optimum is
**no longer a valid lower bound** on \(C_{6.2}\).

(In other words: with the naive (F), v2's "lower bound" can exceed
\(C_{6.2}\). That is the precise sense in which the constraint is
non-rigorous.)

### Issue 2 (FATAL). The repaired (F) has an undefined parameter \(M_2\).

The proof acknowledges Issue 1 and proposes the repaired form
\[
\Omega \;\ge\; \text{partial sum} - 2T_K,
\qquad T_K\ge \sum_{k>K}|a_k^2 - b_k^2| + \sum_{k>K}|2 a_k b_k|.
\]
\(T_K\) is then bounded via Parseval:
\(\sum_k(a_k^2 + b_k^2) \le \int_{-1/4}^{1/4} f^2 \le M_2\) for "any
known upper bound \(M_2\) on \(\int f^2\)".

**There is no a-priori upper bound on \(\int f^2\) in this problem.**
The constraints on \(f\) are only \(f\ge 0\), \(\operatorname{supp}(f)
\subset[-1/4,1/4]\), \(\int f = 1\). \(\int f^2\) is unbounded on
this set: take \(f = (1/\epsilon)\mathbf 1_{[0,\epsilon]}\) (with
\(\epsilon\le 1/4\)), then \(\int f^2 = 1/\epsilon \to \infty\) as
\(\epsilon\to 0\). So the LP has no \(M_2\) to plug in.

The proof leaves \(M_2\) as an external parameter "for any known
upper bound", but offers no construction. Without \(M_2\), the
"rigorous" version of (F) cannot be added to the program.

(One *could* introduce a free LP variable \(M_2\) bounded above by
some other LP-quantity. The natural candidate is \(\Omega\) itself,
via \(\int f^2 \ge \|f*f\|_\infty\) [from \(\|f*f\|_\infty\le \|\hat f
\|_2^2 = \int f^2\) by Plancherel]. But this is a *lower* bound on
\(\int f^2\), not the upper bound (F) needs. There is no LP-side
upper bound on \(\int f^2\). This is a structural obstacle.)

### Issue 3 (FATAL, self-acknowledged). The proof explicitly admits v2's implementation is non-rigorous.

Direct quotes:

> In the present implementation we choose \(K\) large enough that the
> truncation tail is negligible at the working numerical precision and
> simply add the lifted partial-sum inequality without the explicit
> tail slack. **Under that simplification the relaxation may slightly
> violate mathematical rigor by an amount bounded by the truncation
> tail.**

> (The truncation issue is properly handled in v3, which uses a finite
> moment-relaxation that avoids the infinite tail.)

> The v2 (F) test-point inequality is established as above (**modulo
> the tail issue**; see v3 for a fully rigorous treatment).

This is exactly the kind of "approximation that may or may not be a
mathematically valid lower bound" that the verifier instructions
forbid in the strongest possible terms:

> MATHEMATICAL RIGOR IS THE MOST IMPORTANT!!! Do not simply assume
> that you are doing an approximation that is good enough or do not
> propose things that may or may not be a mathematically valid lower
> or upper bound.

A `rigorousproof.md` cannot certify a constraint that the proof
itself disclaims as "may slightly violate mathematical rigor". By the
verifier's own standard, this is INVALID, regardless of whether the
disclaimed amount is "small" in practice.

### Issue 4. Inherited Issue 1 of the v1 verdict is not addressed.

The "Validity of the relaxation" section says:

> The v1 constraints (PSD lift on \(Q\), marginals (B), nonneg (A),
> central window (C), bands) are all satisfied by the rank-one
> \(Q = pp^\top\).

But the v1 band constraints are exactly what my v1 verdict's Issue 1
flagged: the v1 proof's chain
\(\frac12\sum_{j+k\in\{m,m+1\}}Q[j,k] = \int_{B_m}(f*f) \le L\Omega\)
relies on a *false* equality. The LHS is not \(\int_{B_m}(f*f)\) but
\(\int_{B_m}(\bar f*\bar f)\) where \(\bar f\) is the cell-averaged
\(f\); the two are not equal in general, and the inequality
\(\int_{B_m}(\bar f*\bar f) \le L\|f*f\|_\infty\) was never proven
(neither in v0, nor v1, nor v2). The relevant counterexample to the
*equality* — \(f = 4\cdot\mathbf 1_{[-1/8,1/8]}\) on the \(N=1\) grid,
band \(B_1=[-1/2,-1/4]\): LHS \(=1/8\), \(\int_{B_1}(f*f) = 0\) — is
in the v1 verdict.

v2 needed to either fix this or document it as known. It does
neither.

### Issue 5 (cosmetic but worth noting). The "Numerical observation" reveals the SDP lift is degenerate.

The proof itself notes:

> At \(N=8\), \(K=6\), \(n_\text{tests}=41\), MOSEK still returns
> \(\Omega^{(2)} = 1\). The reason is that the SDP lift on \(M\) does
> not constrain the diagonal entries \(M[K+k, K+k]\) from above, so
> the relaxation can artificially inflate \(M[K+k, K+k]\) to make the
> \((M[k,k] - M[K+k, K+k])\cos(2\pi k t_*)\) term arbitrarily negative
> when \(\cos(2\pi k t_*) > 0\), and arbitrarily positive when
> \(\cos(2\pi k t_*) < 0\), balancing out to \(0\) at every \(t_*\).

This is a real structural defect: even if (F) were rigorous, it
wouldn't actually tighten the LP, because the degenerate inflation of
\(M[K+k,K+k]\) makes the constraint trivially satisfied at \(\Omega=1\).
So v2 produces no numerical improvement, *and* the constraint is
non-rigorous, *and* the moment matrix is degenerate. v2 is doubly
broken: it doesn't tighten and it isn't rigorous.

The proposer's plan ("the fix in v3 is to introduce symmetric
three-tensor moments and localizing matrices for \(p_j \ge 0\)") is
the right structural move, but should be paired with the
rigorous-tail treatment of (F).

## 4. What would make this VALID

1. **(F) needs a real tail control.** The cleanest options are:

   * **Replace (F) by per-coefficient constraints.** Instead of
     evaluating the truncated Fourier sum at test points (which
     suffers Gibbs overshoot), constrain the *individual coefficients*
     \(\hat F(k)\) — e.g.\ \(|\hat F(k)| \le \int F = 1\) (trivially)
     and the more useful \(\hat F(k) = \hat f(k)^2\) which is rank-one
     in \((a_k, b_k)\) and matches what (E) already encodes.
   * **Use the de la Vallée-Poussin or Fejér kernel** instead of the
     symmetric partial sum. The Fejér kernel is *nonnegative*, so the
     Fejér-mean is a true lower bound on \(F(t_*)\), avoiding Gibbs
     overshoot. The Fejér mean of order \(K\) is
     \(\sigma_K(t) = \frac{1}{K+1}\sum_{j=0}^{K} S_j(t)\) where
     \(S_j\) are the symmetric partial sums; this is also expressible
     linearly in \((M[k,k], M[K+k,K+k], M[k,K+k])\) and gives
     \(\sigma_K(t_*) \le F(t_*) \le \Omega\) **without** any tail
     correction. This is the right repair.
   * **Or:** use a positive trigonometric polynomial as the test
     function instead of point evaluation. For any nonneg
     trig-polynomial \(P_K\) of degree \(K\),
     \(\int P_K F\ge 0\) and \(\int P_K\,dt\) is computable;
     this gives a *linear* constraint on \(M\) without any tail
     slack.
   
   I strongly recommend the Fejér-kernel or positive-trig-polynomial
   route. Both avoid the Gibbs problem entirely.

2. **Resolve v1's Issue 1 (band constraint).** The right fix is the
   one I suggested in the v1 verdict: bound \(\int_{B_m}(f*f)\) from
   above using the cell-decomposition of the rectangles \(I_j\times
   I_k\) and the *maximum* of the kernel on each rectangle, which
   gives an upper bound directly in terms of \(p_j p_k\) and the
   correct inequality direction. Or equivalently: define \(K_{m,j,k}\)
   as the *area* of the triangular intersection of the support of
   \(\mathbf 1_{I_j}*\mathbf 1_{I_k}\) with the band \(B_m\), and use
   the upper bound
   \(\int_{B_m}(f*f) \le \max_{(s,t)\in I_j\times I_k\,:\,s+t\in B_m} f(s)f(t) \cdot K_{m,j,k}\),
   then majorize the maximum by \(\|f*f\|_\infty\) at the *band scale*
   (not the cell scale). The details need to be worked out and
   written down, but the structure is doable.

3. **Address the SDP degeneracy.** The right fix is the one the
   proof itself suggests for v3: localizing matrices for \(p_j\ge 0\),
   which forces \(M[K+k, K+k]\le M[k,k]\) (or a similar bound) and
   prevents the inflate-cancel trick.

## 5. Feedback to the proposer

The structural direction is right: lifting to a moment matrix indexed
by \((1, a_k, b_k)\) and using the rank-one identification
\(\hat F(k) = \hat f(k)^2\) is exactly the bilinear handle Section 5
of White uses for problem 6.3. The single biggest mistake in v2 is to
test the lifted \(F\) at *point values* via Fourier partial sums
without accounting for Gibbs overshoot. The fix is small (Fejér
kernel instead of partial sum) and bypasses both Issues 1 and 2 in
one move.

The self-acknowledgment of the tail issue in (F) is honest and
appreciated, but per the verifier's standard it cannot be certified
as is. The proof needs to either implement the rigorous tail-bound
form or switch to a Gibbs-free test function.

Suggested next iteration (v3):

1. **Replace (F) by a Fejér-kernel test.** For each test point
   \(t_*\) and each \(K\), define
   \(\sigma_K(t_*) = \frac{1}{K+1}\sum_{j=0}^{K} S_j(t_*)\) and add
   the constraint \(\Omega \ge \sigma_K(t_*)\). This is a finite sum
   of *nonneg-weighted* terms in \(M\), strictly \(\le F(t_*)\le
   \Omega\), with **no tail correction needed**.
2. **Add localizing matrices for \(p_j\ge 0\)** to break the
   diagonal-inflation degeneracy.
3. **Fix v1's band constraint via the cell-rectangle area upper
   bound** (or by switching to Fejér tests on the band integrals as
   well — same idea).
4. **Drop the test-point sweep entirely** in favor of the Fourier-side
   constraint \(\hat F(k) = \hat f(k)^2\), enforced as a rank-one
   relaxation via (E). The rank-one relaxation contains all the
   information of (F) at *every* \(t_*\) simultaneously, without
   overshoot.

Once these are in place, v3 should produce a numerical bound strictly
above 1 (the trivial bound), since the constraints will actually bind.

---

**Final verdict: INVALID.**
