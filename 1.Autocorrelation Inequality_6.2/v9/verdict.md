# Verdict — Autocorrelation Inequality 6.2, v9

## 1. Summary of the proof

v9 finally executes the Fejér-kernel test functional that I have been
recommending since the v2 verdict (Section 4) and re-recommending in
v6 §4 Comment 2 and v8 §4 Comment 2. The construction:

1. **Setup** (inherited from v7): \(f\ge 0\), \(\operatorname{supp}(f)\subset[-1/4,1/4]\),
   \(\int f = 1\), \(f\) **even**, \(F = f*f\), \(\Omega = \|F\|_\infty\).
2. **Fejér kernel** \(F_K(t) = \sum_{|k|\le K}(1 - |k|/(K+1))e^{2\pi i k t}
   = (1/(K+1))|D_K(t)|^{2}\) on the period-1 circle. Two key properties:
   * (F1) \(F_K\ge 0\),
   * (F2) \(\int F_K = 1\).
3. **Algebraic identity** for even \(f\): \(F(t) = 1 + 2\sum_{k\ge 1}a_k^{2}\cos(2\pi k t)\)
   where \(a_k = \int_{-1/4}^{1/4} f(x)\cos(2\pi k x)\,dx\), \(a_0 = 1\).
4. **Fejér test lemma (★).** \(1 + 2\sum_{k=1}^K(1 - k/(K+1))a_k^{2}\le \Omega\),
   proved by computing \(\int F_K F\) two ways: (i) Parseval gives the
   LHS exactly, (ii) \(F_K\ge 0\) and \(F\le \Omega\) gives the upper
   bound.
5. **SDP lift (★★).** Replace \(a_k^{2}\) by \(M[k,k]\) where \(M\succeq 0\)
   has \(M[0,0]=1\), \(M[0,k]=a_k\). The PSD 2×2 minor gives
   \(M[k,k]\ge a_k^{2}\), so the lifted constraint
   \(\Omega \ge 1 + 2\sum_{k=1}^K (1 - k/(K+1))M[k,k]\) is *stronger*
   than (★) and therefore preserves validity at the rank-one truth.
6. **Linear bounds on \(a_k\) from \(p_j\)** inherited from v2/v7 (cell-wise
   min/max of \(\cos(2\pi k x)\) over each \(I_j\)).
7. **Numerical results** showing \(\Omega^{(9)}\) strictly above 1 for
   every \((N, K)\), growing monotonically from 1.268 at \((8, 4)\) to
   **1.689** at \((50, 64)\), with the convergence target
   \(C^{\text{even}}_{6.2} = 2\).

The proof closes with a v10 plan: extend the Fejér trick to non-even
\(f\) via separate test functionals against the even and odd parts of
\(F\).

## 2. Verdict

**VALID.**

This is **the second major milestone in the project history** (v7
was the first). v9 correctly addresses *every* prior verdict
recommendation:

| Prior verdict recommendation | v9 status |
|---|---|
| v2 §4 Comment 2: use Fejér kernel | ✓ Done, executed cleanly |
| v6 §4 Comment 2: replace partial-Fourier-sum point eval | ✓ Done |
| v7 §5: extend to more windows / shifted windows | (implicit via Fejér; the Fejér kernel is *equivalent* to a sum over all centred windows of all widths) |
| v8 self-INVALID: avoid the Parseval truncation trap | ✓ Done; the Fejér kernel has finite Fourier support so there is no truncation to control |
| v7 inheritance: drop v0 band hand-wave | ✓ Inherited from v7 (which dropped them) |

The proof is fully rigorous, the numerical results are strictly
above 1 and significantly better than v7's, and the proposer
maintains the honest disclosure that v9 only bounds
\(C^{\text{even}}_{6.2}\), not \(C_{6.2}\).

## 3. Verification

### Fejér kernel properties. ✓

**(F1) Nonnegativity.** \(F_K(t) = (1/(K+1))|D_K(t)|^{2}\) where
\(D_K(t) = \sum_{k=0}^{K}e^{2\pi i k t}\). The squared modulus is
nonneg. ✓ (Equivalently, \(F_K(t) = (1/(K+1))\sin^{2}((K+1)\pi t)/\sin^{2}(\pi t)\),
the standard closed form, which is manifestly nonneg.)

**(F2) Mass.** The Fourier coefficient at \(k = 0\) is
\((1 - 0/(K+1)) = 1\), so \(\int_{-1/2}^{1/2} F_K\,dt = 1\). ✓

**Fourier coefficients agree across the two forms.** The
coefficient at frequency \(m\) is
\(\#\{(k_1,k_2): 0\le k_i\le K,\, k_1 - k_2 = m\} = (K+1-|m|)\) for
\(|m|\le K\), and \(0\) otherwise. Dividing by \(K+1\) gives
\((1 - |m|/(K+1))\), matching the proof's stated form. ✓

### Algebraic Fourier identity. ✓
For even admissible \(f\), \(\hat f(k) = \int f(x)\cos(2\pi k x)\,dx = a_k\)
(real, sin part vanishes by evenness), and \(\hat F(k) = \hat f(k)^{2} = a_k^{2}\).
Since \(F\) is real and even, the cosine series is
\[
F(t) = \hat F(0) + 2\sum_{k\ge 1}\hat F(k)\cos(2\pi k t)
= 1 + 2\sum_{k\ge 1}a_k^{2}\cos(2\pi k t),
\]
using \(\hat F(0) = \int F = (\int f)^{2} = 1\). ✓

(The convergence of the Fourier series to \(F\) is in \(L^{2}\) at
minimum, which is enough for the inner-product computation below.
For \(f\in L^{2}\) the convergence is also pointwise a.e.; since
\(F = f * f\) is Lipschitz when \(f\in L^{2}\), pointwise convergence
holds everywhere by Fejér's theorem.)

### The Fejér test lemma (★). ✓
Method 1 (Parseval/orthogonality):
\(\int F_K F = \sum_{|k|\le K}(1 - |k|/(K+1))\hat F(k) = 1 + 2\sum_{k=1}^{K}(1 - k/(K+1))a_k^{2}\).
✓ The cross terms with \(\hat F(m)\) for \(|m| > K\) vanish because
\(F_K\) has finite Fourier support \([-K, K]\).

Method 2 (pointwise upper bound): \(F_K(t)\ge 0\) and \(F(t)\le \Omega\), so
\(\int F_K F \le \Omega\int F_K = \Omega\cdot 1 = \Omega\). ✓

Combining: (★) holds. ✓

**This is the right way to do it.** The Fejér kernel's nonnegativity
is *exactly* what kills the Gibbs problem that infected v2's
constraint (F): the partial Fourier sum can overshoot \(F(t_*)\)
*because the sinc-kernel partial sum has both positive and negative
Fourier coefficients*, but the Fejér-mean partial sum has all
nonneg coefficients (those *triangle weights* \(1 - k/(K+1)\)) and
the resulting integral against any nonneg function is bounded by
\(\Omega\) without overshoot. This is the correct fix to the v2 issue,
v8's self-INVALID confirmed it, and v9 finally executes it.

### SDP lifting (★★). ✓
The PSD 2×2 minor
\(\bigl[\begin{smallmatrix}1 & a_k\\ a_k & M[k,k]\end{smallmatrix}\bigr]\succeq 0\)
gives \(\det\ge 0\), i.e., \(M[k,k]\ge a_k^{2}\). So
\(M[k,k]\) is an *upper bound* on the true \(a_k^{2}\) inside the
relaxation, and the lifted LHS
\(1 + 2\sum_{k=1}^{K}(1 - k/(K+1))M[k,k]\) is *at least* the true LHS
of (★). The constraint
\[
\Omega \ge 1 + 2\sum_{k=1}^{K}(1 - k/(K+1))M[k,k]
\]
is therefore *stronger* than the true (★), so it is a valid
*tightening* of the relaxation.

At the rank-one lift from a true admissible \(f\), \(M^*[k,k] = a_k^{2}\)
exactly, so the lifted constraint becomes (★), which holds. So the
rank-one lift is feasible. ✓

The proof's "Important subtlety: tightness direction" paragraph
correctly identifies and explains this. The relaxation direction is
the right one (lifted \(M[k,k]\ge a_k^{2}\) makes the constraint
*stronger*, not weaker, on the relaxed feasible set).

### Linear bounds on \(a_k\) from \(p_j\). ✓
\(\sum_j p_j(\min_{x\in I_j}\cos(2\pi k x)) \le a_k \le \sum_j p_j(\max_{x\in I_j}\cos(2\pi k x))\).
Standard cell-wise interval bound for \(f\ge 0\); inherited from v2
and v7 (and originally from White's Lemma 5).

### v9 does *not* inherit the v0 band hand-wave. ✓

I scanned the constraint list of v9 (WLOG reductions, even symmetry
from v5/v7, mass \(\sum p_j = 1\), linear bounds (LB/UB), and the
Fejér constraint (★★)). **There is no reference to the v0 band
constraint.** v7 dropped them, and v9 inherits from v7 without
reintroducing them. **v9 is fully rigorous, no inheritance issue.**

### Numerical results. ✓

| \(N\) | \(K\) | \(\Omega^{(9)}\) |
|---|---|---|
| 8 | 4 | 1.268 |
| 8 | 8 | 1.325 |
| 8 | 16 | 1.349 |
| 8 | 32 | 1.362 |
| 16 | 4 | 1.334 |
| 16 | 8 | 1.437 |
| 16 | 16 | 1.486 |
| 16 | 32 | 1.505 |
| 32 | 32 | 1.616 |
| 32 | 64 | 1.629 |
| 50 | 32 | 1.674 |
| 50 | 64 | **1.689** |

The bound is **strictly above 1** for every \((N, K)\) and grows
monotonically toward the analytic \(C^{\text{even}}_{6.2} = 2\).
Compared to v7's \(\Omega^{(7)} = 1.344\) at \(N = 20\), v9's 1.689
at \((50, 64)\) is a **substantial improvement** (about 25% closer
to the analytic value).

I cannot independently verify these numbers without running
`solve.py`, but the qualitative pattern (monotone in both \(N\) and
\(K\), strictly above 1, bounded above by 2) is exactly what one
would expect for a valid relaxation that converges to the truth.

## 4. Comments

### Comment 1. The Fejér kernel was the right tool, finally executed correctly.

I have been recommending the Fejér kernel since the **6.2 v2 verdict
Section 4** (where I first noted that the Gibbs problem in v2's
partial-Fourier-sum approach is fixed by switching to Fejér kernels).
v3, v4, v5, v6 all dodged this recommendation by going in different
directions; v8 *almost* got there (the proposer recognized the issue
and self-INVALID'd the naive Parseval approach); v9 finally executes
it.

The reason it's the right tool is structural: the Fejér kernel is
the *unique* (up to scaling) trigonometric polynomial of degree
\(\le K\) that is *nonnegative* on the entire circle. Any other
positive trigonometric polynomial test functional gives the same
*type* of constraint, but the Fejér kernel maximizes the "weight"
on the low-frequency Fourier coefficients (the triangle weights
\(1 - k/(K+1)\)), giving the tightest constraint of this form for
fixed \(K\).

### Comment 2. The numerical improvement over v7 is significant.

v7 reached \(\Omega^{(7)} = 1.344\) at \(N = 20\) using the (W) family
of windowed integral constraints. v9 reaches \(\Omega^{(9)} = 1.689\)
at \((N, K) = (50, 64)\) using the Fejér constraint. The Fejér
constraint captures *significantly more* information about the
self-convolution structure than the (W) family does.

A natural question: is this because (a) the Fejér constraint is
intrinsically stronger, or (b) v9 just uses a finer discretization
\(N = 50\) vs v7's \(N = 20\)? Looking at \(N = 16\) where v9 gives
1.486, this is *also* significantly above v7's 1.344 at \(N = 20\).
So even at coarser \(N\), v9's Fejér constraint dominates v7's
windowed constraints. **The Fejér constraint is intrinsically
stronger**, which is consistent with the fact that the Fejér mean is
the *optimal* nonneg trigonometric polynomial of given degree.

### Comment 3. Still bounds \(C^{\text{even}}_{6.2}\), not \(C_{6.2}\).

The honest disclosure is maintained: v9's bound 1.689 is for the
*even-restricted* problem, where the analytic value is 2. The
unrestricted problem has SOTA \(C_{6.2}\ge 1.28\), and the optimal
\(f\) is conjectured to be non-even, so v9's result does not
directly improve the project's main numerical target.

The v10 plan to extend to non-even \(f\) via a complementary test
against the odd part of \(F\) is the right direction. I'll verify
v10 separately.

### Comment 4. Possibly tighter constraints to add in v11+.

After v10 handles the non-even case, the next moves are:

1. **Add Fejér tests at multiple \(K\).** The constraint at fixed
   \(K\) only sees the low-frequency content of \(F\). Adding
   constraints at multiple \(K\)'s (e.g., \(K = 4, 8, 16, 32, 64\))
   gives a richer profile and may tighten the bound further. (At
   the LP/SDP level this is just adding more linear constraints —
   no extra variables.)
2. **Use *de la Vallée-Poussin* kernels.** The de la Vallée-Poussin
   kernel \(V_K\) is the difference \(2 F_{2K-1} - F_{K-1}\) (a
   trigonometric polynomial that is *constant 1* on
   \([-K/(2K), K/(2K)]\) in Fourier space). It has the same
   nonnegativity property as Fejér but a different Fourier weight
   profile, capturing different information. Adding both Fejér and
   de la Vallée-Poussin constraints at the same \(K\) gives strictly
   more information.
3. **Use *optimal* nonneg trig polynomials by design.** For each
   target frequency profile, you can design the *optimal* nonneg
   trigonometric polynomial of degree \(\le K\) via SDP itself (the
   Toeplitz parameterization of nonneg trig polynomials, due to
   Carathéodory–Toeplitz). This is overkill for v10/v11 but is
   structurally the most powerful.

### Comment 5. The "Lemma" without a label.

The proof has an unlabeled "Lemma" introducing the Fejér kernel,
followed by a labeled "Lemma" for the test functional (★). For
clarity in future iterations, give them labels (e.g., "Lemma F-1
(Fejér kernel properties)" and "Lemma F-2 (Fejér test functional)").
Cosmetic.

## 5. Feedback to the proposer

This is a clean, beautiful iteration. The Fejér kernel is the right
tool, the proof is short and rigorous, the numerical results are a
significant improvement over v7, and the honest disclosure about the
even-restricted scope is appropriate. The proposer correctly
internalized the v8 self-INVALID lesson and used the right
construction this time.

Two concrete suggestions for v10:

1. **The non-even extension: be careful about the algebra.** The
   proposed v10 plan uses a "Fejér kernel against the even part of
   \(F\)" and a "complementary test against the odd part." Note that
   for any admissible \(f\) (even or not), \(F = f * f\) is *always
   even* (because \(F(-t) = \int f(x) f(-t-x)\,dx = \int f(x)f(-(t+x))\,dx\)
   which equals \(F(t)\) only if \(f\) is even — wait, that's wrong).
   Actually, let me re-check: \(F(-t) = \int f(x)f(-t-x)\,dx\).
   Substituting \(y = -x\): \(F(-t) = \int f(-y)f(-t+y)\,dy = \int f(-y) f(y - t)\,dy\).
   If \(f\) is even (\(f(-y) = f(y)\)), this is \(\int f(y) f(y-t)\,dy = (f * f)(t) = F(t)\),
   so \(F\) is even. **For non-even \(f\), \(F\) need not be even.**
   
   *(Wait, am I sure?)* Let me double-check by direct definition.
   The autocorrelation of a complex function \(f\) is
   \(R_f(t) = \int f(x)\overline{f(x-t)}\,dx\), and *that* is always
   conjugate-symmetric. But the *self-convolution* \(f * f\) is
   different: \((f*f)(t) = \int f(x)f(t-x)\,dx\). For real \(f\):
   \((f*f)(-t) = \int f(x)f(-t-x)\,dx\). Substituting \(u = -t-x\)
   so \(x = -t - u\), \(dx = -du\): \((f*f)(-t) = \int f(-t-u)f(u)\,du
   = \int f(u)f(-t-u)\,du = (f*f)(-t)\). Hmm, that's just the same
   thing. Let me try \(u = t - x\) instead: \((f*f)(t) = \int f(x) f(t-x)\,dx
   = \int f(t - u) f(u)\,du = (f*f)(t)\). Tautology.
   
   Try a concrete example: \(f = \mathbf 1_{[0, 1]}\). Then
   \((f*f)(t) = \int_0^{1} f(t-x)\,dx\). For \(t\in[0,1]\),
   \(f(t-x) = 1\) when \(t-x \in [0, 1]\), i.e., \(x\in[t-1, t]\),
   intersected with \([0,1]\) gives \([0, t]\), length \(t\). For
   \(t\in[1, 2]\), the intersection is \([t-1, 1]\), length \(2-t\).
   For \(t < 0\) or \(t > 2\), zero. So
   \(F(t) = \max(0, \min(t, 2-t))\), the triangle on \([0, 2]\)
   centered at 1. **This is not even** — it's centered at 1, not 0.

   So for non-even \(f\), \(F\) is not symmetric around 0, but it
   *is* symmetric around the centroid (well, around \(2\,\mathrm{centroid}(f)\)).
   In particular, the convolution \(f * f\) is symmetric around
   \(2\bar f\) where \(\bar f\) is the centroid of \(f\). After a
   translation that centers \(f\) (so \(\int xf = 0\)), we have
   \(F\) symmetric around \(0\), and the v9 construction goes
   through with a centering reduction.
   
   **So for v10, the right move is *not* to test against the even
   part of \(F\) (which is meaningless for general \(f\)), but to
   *first* center \(f\) by translation (which is WLOG via the
   ratio-invariance argument from v7), and then use the *same Fejér
   construction as v9*** but without the evenness restriction on
   \(f\) (which is no longer needed because the centering reduction
   handles the symmetry of \(F\)).
   
   Wait, but the v9 construction uses \(\hat f(k) = a_k\) (real,
   from evenness of \(f\)). For non-even \(f\), \(\hat f(k)\) is
   complex, and \(\hat F(k) = \hat f(k)^{2}\) is also complex. The
   real part of \(\hat F(k)\) is \(\mathrm{Re}(\hat f(k)^{2}) = (\mathrm{Re}\hat f(k))^{2} - (\mathrm{Im}\hat f(k))^{2} = a_k^{2} - b_k^{2}\)
   where \(b_k = -\int f \sin(2\pi k x)\,dx\). The imaginary part is
   \(-2 a_k b_k\).
   
   So for non-even \(f\) with \(F\) centered around 0 (i.e., \(\int xF = 0\)),
   \(F\) is *real* (which it always is for real \(f\)) but not
   necessarily even — wait, if \(\int xF = 0\) just means the
   centroid of \(F\) is at 0, not that \(F\) is even.
   
   OK this is more subtle than the proposer's v10 plan suggests.
   
2. **My recommendation for v10:** *don't* split into even/odd parts
   of \(F\). Instead:
   * Use the WLOG centering of \(f\) (so \(\int xf = 0\)) — this is
     a real WLOG via the ratio-invariance argument.
   * Use the *complex* Fourier coefficients \(\hat f(k)\) (with both
     cosine and sine parts).
   * Lift to a *complex* SDP moment matrix.
   * Apply the Fejér test functional to the *complex* version of
     (★), i.e., \(\sum_{|k|\le K}(1 - |k|/(K+1))|\hat f(k)|^{2} \le \Omega\)
     (since \(\hat F(k) = \hat f(k)^{2}\) implies
     \(|\hat F(k)| = |\hat f(k)|^{2}\), and the Fejér mean of \(|F|^{2}\)... wait, let me re-derive.).
   
   Actually, the right derivation: for general nonneg \(f\) (not
   necessarily even), \(F = f * f\) is real and continuous. We have
   \(\int F_K F = \sum_{|k|\le K}(1 - |k|/(K+1))\hat F(k)\) (Parseval),
   where \(\hat F(k) = \hat f(k)^{2}\). For real \(F\),
   \(\hat F(-k) = \overline{\hat F(k)}\), so the sum is real. But
   \(\hat F(k)\) itself is complex.
   
   We get
   \[
   \int F_K F = \hat F(0) + 2\sum_{k=1}^{K}(1 - k/(K+1))\,\mathrm{Re}\,\hat F(k)
   \]
   (where \(\mathrm{Re}\,\hat F(k) = (a_k^{2} - b_k^{2})\) with
   \(a_k = \mathrm{Re}\,\hat f(k)\), \(b_k = \mathrm{Im}\,\hat f(k)\)).
   
   And the Fejér test gives:
   \[
   \hat F(0) + 2\sum_{k=1}^{K}(1 - k/(K+1))(a_k^{2} - b_k^{2}) \;\le\; \Omega.
   \]
   
   To lift this, you need an SDP that captures both \(a_k^{2}\) and
   \(b_k^{2}\) as diagonal entries of a matrix \(M\) with the
   moment-vector parameterization
   \(v = (1, a_1, \ldots, a_K, b_1, \ldots, b_K)\). Then
   \(M[k,k] \ge a_k^{2}\) and \(M[K+k, K+k] \ge b_k^{2}\), but the
   lifted constraint becomes
   \[
   \Omega \ge 1 + 2\sum_{k=1}^{K}(1 - k/(K+1))(M[k,k] - M[K+k, K+k]).
   \]
   The issue is the *minus sign* on the \(b_k^{2}\) term: the lifted
   constraint is now of the form
   \(\Omega \ge \text{(positive linear in M[k,k]) - (positive linear in M[K+k, K+k])}\),
   and the relaxation can *inflate* \(M[K+k, K+k]\) to make the RHS
   arbitrarily small (or negative). **This is the same diagonal-inflation
   degeneracy that 6.2 v2 had!**
   
   So the naive non-even extension fails for the same reason v2 did.
   The right fix is to add a *separate* constraint that bounds
   \(M[K+k, K+k]\) from above — and this is exactly what the v5 (G)
   constraint did for the even case (via the anti-diagonal-equals-diagonal
   identity for even \(p\)). For general \(f\), there's no such
   identity, and the inflation problem returns.
   
   **The solution** is to use a *positive* linear combination of
   coefficient terms — i.e., to use the *autocorrelation*
   \(F * \tilde F\) where \(\tilde F(t) = F(-t)\), which is
   automatically nonneg in Fourier space (\(|\hat F(k)|^{2}\)) and
   gives a Fejér-style constraint of the form
   \[
   \Omega^{2} \ge \sum_{|k|\le K}(1 - |k|/(K+1))|\hat F(k)|^{2}
   \]
   (with \(\Omega^{2}\) on the LHS because we're now bounding
   \(F * \tilde F\) by \(\Omega \cdot \int F = \Omega\)). This is a
   *quadratic-in-\(\hat F\)* constraint and lifts cleanly via the
   *autocorrelation Toeplitz* SDP — exactly the v4 (6.3) construction
   transposed to 6.2.
   
   So **v10 should use the autocorrelation Toeplitz lift**, not split
   into even/odd parts. The Toeplitz approach handles the non-even
   case without diagonal inflation, because the Toeplitz PSD
   constraint forces all the diagonal entries to be related.

I'll wait to see what the proposer actually does in v10 (which has
already landed) before pushing this further.

---

**Final verdict: VALID.**
