# Verdict — Autocorrelation Inequality 6.2, v54

## 1. Summary

v54 adds a second localizer
\[
h_2(x) \;:=\; \cos(2\pi x) - \cos^{2}(2\pi x) \;=\; \cos(2\pi x)\bigl(1 - \cos(2\pi x)\bigr)
\]
to v51's setup. \(h_2\) is nonneg on \([-1/4, 1/4]\) (same support
region as v51's \(h_1 = \cos(2\pi x)\)), but algebraically
different. The corresponding localized Bochner Toeplitz
\(T_{\nu_2}\succeq 0\) is added, with
\[
\hat\nu_2(k) = -\tfrac{1}{2}\hat f(k) + \tfrac{1}{2}(\hat f(k-1) + \hat f(k+1)) - \tfrac{1}{4}(\hat f(k-2) + \hat f(k+2)),
\]
requiring Fourier moments up to index \(K + 2\).

v54 = v51 + \(T_{\nu_2}\succeq 0\).

## 2. Verdict

**VALID** (at finite truncation; in the infinite limit, the
constraint is implied by v51's \(h_1\) constraint via SOS, but at
finite Toeplitz size it is a genuinely new PSD cut).

## 3. Verification

### (1) Nonnegativity on \([-1/4, 1/4]\). ✓

With \(c = \cos(2\pi x)\):
* On \([-1/4, 1/4]\), \(c\in[0, 1]\), so \(c(1 - c)\in[0, 1/4]\)
  — nonneg.
* On \((1/4, 3/4)\), \(c\in[-1, 0)\), so \(c(1-c) < 0\) (since \(c < 0\), \(1 - c > 0\)) — negative.
* Maximum \(c(1 - c) = 1/4\) attained at \(c = 1/2\), i.e.,
  \(x = \pm 1/6\).

So \(h_2\ge 0\) on \([-1/4, 1/4]\), matching the support of \(f\). ✓

### (2) Fourier expansion. ✓

\[
\cos^{2}(2\pi x) = \tfrac{1 + \cos(4\pi x)}{2}
\]
\[
h_2 = \cos(2\pi x) - \tfrac{1 + \cos(4\pi x)}{2}
= -\tfrac{1}{2} + \cos(2\pi x) - \tfrac{1}{2}\cos(4\pi x).
\]
Fourier coefficients:
* \(\hat h_2(0) = -1/2\) ✓
* \(\hat h_2(\pm 1) = 1/2\) ✓
* \(\hat h_2(\pm 2) = -1/4\) ✓
* All others zero. ✓

### (3) Convolutional computation of \(\hat\nu_2\). ✓

\[
\hat\nu_2(k) = \widehat{h_2 f}(k) = \sum_m \hat h_2(m)\,\hat f(k - m)
\]
\[
= -\tfrac{1}{2}\hat f(k) + \tfrac{1}{2}(\hat f(k-1) + \hat f(k+1)) - \tfrac{1}{4}(\hat f(k-2) + \hat f(k+2)).
\]
✓. The real and imaginary parts with the sign-of-index factors
are correctly stated.

### (4) Validity of the relaxation.

For any admissible \(f\ge 0\) supported in \([-1/4, 1/4]\),
\(\nu_2 = h_2 f\) is a positive measure on \(\mathbb T\) (product of
two nonneg functions), and by Bochner, \(T_{\nu_2}\succeq 0\) at
any Toeplitz size. At the rank-one truth, the Fourier moments give
\(T_{\nu_2}\) explicitly, and the constraint is satisfied.

\(\Omega^{(54)}\ge \Omega^{(51)}\) provably. ✓

### (5) The SOS subtlety: finite vs infinite truncation.

In the infinite-dimensional Schmüdgen hierarchy, the constraint
\(T_{\nu_2}\succeq 0\) is **implied** by \(T_{\nu_1}\succeq 0\)
(v51's constraint), via an SOS rewriting:

\[
h_2 = h_1\cdot (1 - h_1) = h_1\cdot 2\sin^{2}(\pi x),
\]
and \(\sin^{2}(\pi x)\) is SOS (it's \((\sin\pi x)^{2}\)). So for
any test polynomial \(\phi\),
\[
\int h_2 f |\phi|^{2}\,dx = \int h_1 f \cdot |\sqrt 2 \sin(\pi x)\cdot \phi|^{2}\,dx,
\]
which is a "\(h_1 f\) inner product of a different polynomial." If
the infinite Toeplitz \(T_{\nu_1}^{\infty}\succeq 0\), so is
\(T_{\nu_2}^{\infty}\) — they differ by an SOS multiplier.

**But at finite truncation size \(K + 1\)**, the SOS rewriting
requires test polynomials of degree up to \(K + 1\) (since
\(\sin(\pi x)\cdot \phi\) has one more frequency than \(\phi\)).
So the finite-size \(T_{\nu_2}\) of size \(K + 1\) is NOT
automatically PSD given \(T_{\nu_1}\) of size \(K + 1\) — it would
require \(T_{\nu_1}\) of size \(K + 2\) or larger (via the SOS
expansion).

So at finite truncation, v54's \(T_{\nu_2}\succeq 0\) is a
**genuinely new** constraint at size \(K + 1\).

The proposer correctly notes the "difference of PSD matrices"
observation but doesn't explicitly state the SOS relationship.
Adding that would clarify the claim.

## 4. Comments

### Comment 1. v53 vs v54: different families.

* v53 uses \(h_3 = \cos^{3}(2\pi x)\) — odd powers of \(h_1\).
* v54 uses \(h_2 = \cos(2\pi x)(1 - \cos(2\pi x))\) — \(h_1\)
  times an SOS \((1 - h_1)\).

Both are valid Schmüdgen localizers in the preordering generated
by \(h_1\). In the infinite limit, all such constraints are
asymptotically redundant with \(T_{\nu_1}^{\infty}\succeq 0\), but
at finite truncation each provides a distinct PSD cut.

**v53 and v54 should be combined** into a single iteration that
adds both \(h_3\) and \(h_2\) simultaneously. There's no reason
to split them.

### Comment 2. Expected gain.

At finite truncation, v54's \(T_{\nu_2}\succeq 0\) at size \(K + 1\)
should give a small but nonzero improvement over v51. Order of
magnitude: \(\sim 10^{-4}\), based on the v45/v46 sequence's
pattern of small gains per additional PSD block.

Combined with v53 and v51's original \(h_1\), the three
Schmüdgen localizers may give \(\sim 3\times 10^{-4}\) improvement.
Not game-changing, but real.

### Comment 3. Compute scaling.

\(T_{\nu_2}\) is \((K+1)\times(K+1)\) — same size as \(T_{\nu_1}\).
Adding \(T_{\nu_2}\succeq 0\) adds \(\sim (K+1)^{2}\) entries to
the SDP, comparable to the v51 Bochner block. Total cost
increases by \(\sim 2\times\) over v51.

### Comment 4. Inheritance is clean.

v54 inherits from v51 and extends aux Fourier variables to
\(K + 2\). No new bugs.

## 5. Feedback to the proposer

For v55:

1. **Combine v53 and v54** into a single iteration. Both add
   Schmüdgen localizers and should be run together.
2. **Run v51/v53/v54 combined** at sanity and production
   parameters. Tabulate \(\Omega^{(53+54)} - \Omega^{(51)}\).
3. **Don't explore further localizers** until v51 + v53 + v54 is
   shown to give a meaningful improvement. The Schmüdgen
   hierarchy has diminishing returns quickly.
4. **Focus on what's binding**: the v44 ceiling at \(K = 96\) is
   \(\sim 1.294\). If v51 + v53 + v54 reaches \(\ge 1.295\), the
   project has a new milestone. If not, the focus should return
   to pushing \(K\) or combining with v50's Gram constraint.

---

**Final verdict: VALID** (Schmüdgen-style localizer
\(h_2 = \cos(1 - \cos)\) added to v51; math correctly derived
(triple-identity, Fourier expansion, convolutional moment
computation, rank-one validity). At finite truncation, the
constraint is a genuinely new PSD cut; in the infinite limit, it
is implied by v51's \(h_1\) constraint via an SOS rewriting
involving \(\sin^{2}(\pi x) = (1 - h_1)/2\). Expected gain is
modest (\(\sim 10^{-4}\)). Best bound on \(C_{6.2}\) remains v44's
\(1.2913\)).

---

## 6. Second-pass critical addenda (independent re-review)

A fresh, skeptical re-reading of `rigorousproof.md` confirms the
VALID verdict but surfaces four presentational / pedagogical
issues that should be fixed before v55:

**(A) The SOS-on-the-torus must use integer-frequency polynomials.**
Section 5 of this verdict correctly identifies the SOS implication
but writes \(1 - h_1 = 2\sin^{2}(\pi x)\). Strictly, \(\sin(\pi x)\)
is *antiperiodic* on \(\mathbb T = \mathbb R/\mathbb Z\) and is **not** a
trigonometric polynomial (half-integer frequency). The correct
SOS on \(\mathbb T\) is
\[
1 - \cos(2\pi x) \;=\; \tfrac{1}{2}|1 - e^{2\pi i x}|^{2},
\]
i.e. \(1 - h_1 = |p(z)|^{2}/2\) with \(z=e^{2\pi i x}\) and
\(p(z) = 1 - z\) of degree \(1\). The "test polynomial of degree
up to \(K+1\)" claim then becomes precise: multiplying
\(\phi(z) = \sum_{|k|\le K} c_k z^k\) by \(p(z)/\sqrt{2}\) gives a
Laurent polynomial of degree \(K+1\). This is the shift that
makes \(T_{\nu_2}^{K+1}\) genuinely new given only
\(T_{\nu_1}^{K+1}\), and implied given \(T_{\nu_1}^{K+2}\).

**(B) The "not a linear combination of nonneg Bochner summands"
argument in rigorousproof.md is weak.** Writing
\(T_{\nu_2}=T_{\nu_1}-T_{\nu_{c^2}}\) and noting "difference of
PSDs ≠ PSD" is correct but does not establish that the
constraint is *new* relative to the Schmüdgen preordering — it
only establishes that the constraint is not trivially derivable
by adding two PSD constraints. The actually relevant fact is
point (A): at finite \(K\), \(T_{\nu_2}^{K+1}\) is *not* implied
by \(T_{\nu_1}^{K+1}\) alone because the SOS lift raises the
degree. Section "Why $h_2$ is not a linear combination …" in
`rigorousproof.md` should be rewritten with this framing.

**(C) The ambient reduction "\(\mathrm{supp}\,f\subset[-1/4,1/4]\)"
is assumed silently.** The original problem is on \(\mathbb R\),
not on \(\mathbb T\), and with unconstrained support. The
standard reduction is: replace \(f\) by \(\tilde f = f\cdot\mathbf{1}_{[-1/4,1/4]}\);
then \(N(\tilde f)\le N(f)\) (convolution pointwise decreases)
and \(D(\tilde f) = D(f)\), so \(R(\tilde f)\le R(f)\), hence
\(\inf R = \inf_{\tilde f} R\) over \([-1/4,1/4]\)-supported
functions. Since \(\tilde f\) is supported in
\([-1/4,1/4]\subset(-1/2,1/2)\), its autocorrelation on the
circle agrees with that on the line on \(|t|\le 1/2\), and Fourier
analysis on \(\mathbb T\) applies. This reduction is the
backbone of both v51 and v54 but is nowhere stated in
`rigorousproof.md`. Make it explicit.

**(D) "Not in the span of existing ones" is not the same as
"improves the SDP."** The Status paragraph claims the new
constraint "is not in the span of the existing ones." Linear
span is irrelevant — what matters is whether the PSD cone
\(\{T_{\nu_2}\succeq 0\}\) is redundant given the other PSD
cones. That is a cone/face question, not a linear-span
question. Reword.

**Bottom line:** The mathematics of v54 is correct and the
relaxation is genuinely tightening at finite \(K\). The
presentation in `rigorousproof.md` is sloppy around (A)–(D).
VALID stands.

