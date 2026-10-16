# Verdict — Autocorrelation Inequality 6.2, v52

## 1. Summary

v52 extends v51's localized Bochner Toeplitz from size \((K+1)\)
to size \((2K+1)\), analogously to how v46 extended v45's plain
Bochner. The construction:

* Extend the \(a, b\) variables to indices \(0, 1, \dots, 2K+1\)
  (one past the v46 extension, because \(\hat\nu(k) = (\hat f(k-1) + \hat f(k+1))/2\)
  needs \(\hat f\) one index higher than the target \(\hat\nu\) index).
* Plain Bochner Toeplitz at size \((2K+2)\) on the extended
  \((a, b)\) (from v46).
* Localized Bochner \(T_\nu^{(2K)}\) of size \((2K+1)\) via
  \(\hat\nu(k) = (\hat f(k-1) + \hat f(k+1))/2\) for \(k = 0, \dots, 2K\).

No numerical results yet. The proposer references "the huge jump
from v45 (1.22) to v51 (1.27) at small \(K\)," suggesting v51
sanity-check numbers may have come in (but I haven't seen them
directly).

## 2. Verdict

**VALID** — straightforward extension of v51's localized Bochner
to a larger Toeplitz size; same math at a larger matrix
dimension.

## 3. Verification

Validity inherits from v51 (verified in v51's verdict). The
larger Toeplitz \(T_\nu^{(2K)}\succeq 0\) follows from Bochner
applied to the positive measure \(\nu = h f\) (with
\(h = \cos(2\pi x)\ge 0\) on \([-1/4, 1/4]\)) at size \((2K+1)\).

The proof is the same \(z^{*} T_\nu^{(2K)} z = \int hf|\phi|^{2}\ge 0\)
one-liner; only the matrix dimension grows.

\(\Omega^{(52)}\ge \Omega^{(51)}\) provably (same relaxation
argument).

## 4. Comments

### Comment 1. Analogous to the v45 → v46 extension.

* v45 (Bochner at size \(K+1\)) → v46 (Bochner at size \(2K+1\)) ≈
  small gain (\(\sim 10^{-5}\)).
* v51 (localized Bochner at size \(K+1\)) → v52 (localized Bochner
  at size \(2K+1\)) ≈ ?

If v51's localized Bochner is more binding than v45's plain
Bochner (as the proposer's "1.22 → 1.27" remark suggests), then
v52's extension may give a larger relative improvement than
v46's extension. But this is speculative until numerical results
come in.

### Comment 2. The "1.22 → 1.27" claim needs context.

The proposer's passing remark that v51 gave "1.27" at small \(K\)
while v45 gave "1.22" at small \(K\) is a **substantial jump**
(\(+0.05\)) if true. But:

* v45 at small \(K\) (sanity table: \(K = 8, 16, 32\)) gave bounds
  \(\sim 1.08\)–\(1.22\), all consistent with v44. A jump to 1.27
  at the same parameters would be a major result.
* The proposer hasn't shown v51's sanity table yet, so the "1.27"
  claim is unverified.
* If true, v51 would be the project's new best framework, and
  v52's extension is worth running.
* If the "1.22 → 1.27" was a typo or different parameter
  setting, the claim should be retracted.

The proposer should produce the v51 sanity table to confirm.

### Comment 3. Compute scaling.

The \(T_\nu^{(2K)}\) block is \((2K+1)\times(2K+1)\), so \(\sim 4\times\)
the v51 block. Combined with v46's extended plain Bochner
(\((2K+2)\times(2K+2)\)), v52's total SDP block size is
substantial but still tractable at moderate \(K\).

### Comment 4. Inheritance is clean.

v52 inherits from v51 + v46. The extended aux variables are
consistently defined.

## 5. Feedback to the proposer

For v53:

1. **Confirm v51's sanity numbers** before running v52. If v51
   really jumps from v45's 1.22 to 1.27 at small \(K\), that's
   the main result and v52 is incremental.
2. **Run v51 and v52 production sweep**: \((1000, 64), (1000, 96), (3000, 96)\).
   If v51 alone gives \(\ge 1.30\), that's a major milestone.
3. **Combine v51 with v44's Fejér–Riesz** — they should be
   orthogonal (v44 tightens the constraint side, v51 tightens the
   \(f\)-moment side).
4. **Tabulate the improvements** explicitly so cross-iteration
   comparisons are unambiguous.

---

**Final verdict: VALID** (rigorous extension of v51's localized
Bochner to Toeplitz size \(2K+1\); same math as v51 at a larger
matrix dimension. \(\Omega^{(52)}\ge \Omega^{(51)}\) provably; no
numerical results yet. The proposer's claimed "1.22 → 1.27" jump
from v45 to v51 at small \(K\) needs confirmation via a sanity
table. Best bound on \(C_{6.2}\) remains v44's \(1.2913\)).
