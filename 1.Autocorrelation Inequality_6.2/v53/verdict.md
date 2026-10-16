# Verdict — Autocorrelation Inequality 6.2, v53

## 1. Summary

v53 adds a second Putinar localizer to v51's construction:
\[
h_3(x) \;:=\; \cos^{3}(2\pi x) \;=\; \tfrac{3}{4}\cos(2\pi x) + \tfrac{1}{4}\cos(6\pi x),
\]
which is nonneg on \([-1/4, 1/4]\) (since \(\cos^{3}\) preserves
the sign of \(\cos\), which is nonneg there).

The localized Bochner \(T_{\nu_3}\succeq 0\) is added, with Fourier
moments
\(\hat\nu_3(k) = \tfrac{3}{8}(\hat f(k-1) + \hat f(k+1)) + \tfrac{1}{8}(\hat f(k-3) + \hat f(k+3))\),
requiring Fourier variables up to index \(K + 3\).

v53 = v51 + \(T_{\nu_3}\succeq 0\).

## 2. Verdict

**VALID** — standard higher-level Schmüdgen localizer; soundly
derived.

## 3. Verification

### (1) The \(\cos^{3}\) identity. ✓

Standard triple-angle formula:
\[
\cos^{3}\theta = \cos\theta\cdot\cos^{2}\theta
= \cos\theta\cdot\tfrac{1 + \cos 2\theta}{2}
= \tfrac{\cos\theta}{2} + \tfrac{\cos\theta\cdot\cos 2\theta}{2}
\]
\[
= \tfrac{\cos\theta}{2} + \tfrac{\cos\theta + \cos 3\theta}{4}
= \tfrac{3\cos\theta}{4} + \tfrac{\cos 3\theta}{4}. \;\checkmark
\]

With \(\theta = 2\pi x\): \(\cos^{3}(2\pi x) = \tfrac{3}{4}\cos(2\pi x) + \tfrac{1}{4}\cos(6\pi x)\). ✓

### (2) Nonnegativity on \([-1/4, 1/4]\). ✓

On \([-1/4, 1/4]\), \(2\pi x\in[-\pi/2, \pi/2]\), so \(\cos(2\pi x)\ge 0\).
Cubing preserves the sign: \(\cos^{3}(2\pi x)\ge 0\) on \([-1/4, 1/4]\).
(Outside, \(\cos\) is negative, so \(\cos^{3}\) is also negative.) ✓

**Important observation:** \(\cos(6\pi x)\) *alone* is **not** a
valid localizer — on \([-1/4, 1/4]\), \(6\pi x\in[-3\pi/2, 3\pi/2]\),
and \(\cos\) is not nonneg on this range. But the specific linear
combination \(\tfrac{3}{4}\cos(2\pi x) + \tfrac{1}{4}\cos(6\pi x) = \cos^{3}(2\pi x)\)
is nonneg. So v53's constraint is not trivially implied by plain
v51.

### (3) Fourier moments of \(\nu_3 = h_3 f\). ✓

\[
\hat\nu_3(k) = \int f(x)\Bigl(\tfrac{3}{4}\cos(2\pi x) + \tfrac{1}{4}\cos(6\pi x)\Bigr)e^{-2\pi i k x}\,dx
\]
\[
= \tfrac{3}{4}\cdot \tfrac{1}{2}(\hat f(k-1) + \hat f(k+1))
+ \tfrac{1}{4}\cdot \tfrac{1}{2}(\hat f(k-3) + \hat f(k+3))
\]
\[
= \tfrac{3}{8}(\hat f(k-1) + \hat f(k+1)) + \tfrac{1}{8}(\hat f(k-3) + \hat f(k+3)). \;\checkmark
\]

Note: the real-imaginary decomposition with the \(\text{sgn}\)
factors handles the negative-index cases consistently, as in v51.

### (4) Validity of the relaxation.

For any admissible \(f\ge 0\) supported in \([-1/4, 1/4]\),
\(\nu_3 = h_3 f\) is a positive measure on \(\mathbb T\) (product of
nonneg on the support, zero off the support). By Bochner,
\(T_{\nu_3}\succeq 0\) at any size. Rank-one truths from admissible
\(f\) satisfy \(T_{\nu_3}\succeq 0\), so the relaxation is valid.

\(\Omega^{(53)}\ge \Omega^{(51)}\) provably. ✓

### (5) Why v53 is not redundant with v51 or with higher-\(K\) Bochner.

v51's localizer \(h_1 = \cos(2\pi x)\) gives a Toeplitz constraint
involving \(\hat f(k\pm 1)\). v53's \(h_3 = \cos^{3}(2\pi x)\)
gives a Toeplitz constraint involving both \(\hat f(k\pm 1)\) and
\(\hat f(k\pm 3)\). The two constraints are **independent** (they
are not linear combinations of each other), and the Schmüdgen
localized-moment hierarchy treats them as distinct tightening
layers.

The proposer correctly notes:
* \(h_1\) and \(h_3\) are not proportional.
* \(h_3\) couples frequencies \(k\pm 1\) AND \(k\pm 3\).
* v53 requires Fourier moments up to index \(K + 3\) (one more
  than v51's \(K + 1\)).

## 4. Comments

### Comment 1. Schmüdgen vs Putinar.

The Putinar framework uses localizers \(\{g_i\}\) that collectively
bound the support. The **Schmüdgen** framework extends this by
allowing *products* of localizers. \(h_3 = h_1^{3}\) is an example
of the Schmüdgen extension: higher powers of the same localizer.

For the canonical Putinar setup with localizer
\(g = \cos(2\pi x)\), the full Schmüdgen module includes
\(g, g^{2}, g^{3}, \dots\). But:
* \(g^{2} = \cos^{2}\) is nonneg *everywhere*, so its localized
  Bochner is just the plain (non-localized) Bochner at shifted
  frequencies — no new information beyond v45.
* \(g^{3} = \cos^{3}\) has the same sign as \(g\) (nonneg on
  \([-1/4, 1/4]\), negative elsewhere), so it's a valid localizer
  *for the support restriction*. This gives genuinely new
  information.

v53 uses this directly.

### Comment 2. Higher odd powers.

The same argument applies to \(g^{5} = \cos^{5}\), \(g^{7}\), etc.
Each is a degree-\(k\) trig polynomial with \(\cos\) and
\(\cos 3, \cos 5, \dots, \cos k\) terms, and each gives an
independent localized Bochner constraint requiring Fourier
moments up to \(k + K\).

v54+ could add \(g^{5}\) as another layer, though the gains
from successive odd powers probably diminish.

### Comment 3. Compute cost.

The Toeplitz \(T_{\nu_3}\) is size \((K+1)\times(K+1)\) — same size
as v51's \(T_{\nu_1}\). The only extra cost is the two additional
Fourier variables \(a_{K+2}, a_{K+3}, b_{K+2}, b_{K+3}\), with
their cell bounds. Overall cost is \(\sim 10\%\) more than v51.

### Comment 4. Expected gain.

Small but nonzero. At small \(K\), the additional localizer gives
information about frequencies \(k\pm 3\) that wasn't in v51. At
larger \(K\), the information is already implicitly encoded by
v46's higher-degree plain Bochner + v51's \(h_1\) localizer, so
the marginal gain likely shrinks.

Realistic estimate: \(\sim 10^{-4}\) at moderate parameters.

### Comment 5. Inheritance is clean.

v53 inherits from v51 and extends the auxiliary variables by two
indices. No new bugs.

## 5. Feedback to the proposer

For v54:

1. **Run v53 at sanity parameters** and verify \(\Omega^{(53)}\ge \Omega^{(51)}\).
2. **Run v53 at production parameters** and compare to v51.
   If v51's jump to \(\sim 1.27\) (from the proposer's earlier
   remark) is real, and v53 adds another \(\sim 10^{-3}\)–\(10^{-4}\),
   then v51 + v53 approaches 1.30.
3. **Before running v54 (higher localizer powers)**, confirm v53
   is not vacuous. If it is, the localizer hierarchy is
   saturating.
4. **Consolidate**: v51 (one localizer) + v53 (two localizers)
   should be a single iteration, not two. Going forward,
   combine related localizer additions into one document.

---

**Final verdict: VALID** (higher-power Schmüdgen localizer
\(h_3 = \cos^{3}(2\pi x)\) added to v51; the math is standard
(triple-angle formula, Fourier moments, Bochner at the rank-1
truth). \(\Omega^{(53)}\ge \Omega^{(51)}\) provably. Expected
gain is small — maybe \(10^{-4}\) — because most of the
support-restriction information is already captured by v51 and
v46's plain Bochner at higher \(K\)).
