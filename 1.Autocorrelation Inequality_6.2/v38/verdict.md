# Verdict — Autocorrelation Inequality 6.2, v38

## 1. Summary

v38 is a **project-final summary** rather than a new iteration. No
new mathematical content, no new constraints, no new numerical
runs. The proposer declares the project complete with the
v35 bound:

\[
\boxed{C_{6.2} \;\ge\; 1.2851 \quad \text{(rigorous, v35 at $(N, K) = (3000, 64)$).}}
\]

The summary catalogs the project's history, the v30 construction,
and lists open work directions for future iteration.

## 2. Verdict

**VALID** as a project summary. No new claims to verify; the
construction and validity were established in v30 and the numerical
result in v35.

## 3. Verification

* The claimed bound \(1.2851\) is the v35 result, which I verified
  in `v35/verdict.md`.
* The v30 construction summary in §"Key construction" matches the
  actual v30 LP I verified in `v30/verdict.md`.
* The validity arguments cited (Hölder chain, Fejér kernel
  positivity, no even-restriction) are correct.
* The project history table is consistent with the verdict files
  I have written, *with two minor inaccuracies* — see Comment 1.

## 4. Comments

### Comment 1. Two minor errors in the project history table.

* Row "v15–v23 | combinations | 1.12 | none | VALID (no improvement)":
  The row inappropriately lumps in v22, which was self-INVALID
  (proposer caught the midpoint-quadrature error). The next row
  corrects this for v22 but the v15–v23 row should be
  "v15–v23 (excluding v22)".
* Row "v36–v38 | refinements": this is misleading because v36 was
  a vacuous constraint addition (Bochner Toeplitz) and v37/v38 are
  in-progress / summary, not refinements. None contributed any
  numerical improvement.

These are cosmetic issues; the bottom-line bound is correctly
attributed to v35.

### Comment 2. The proposer declared completion without running the recommended \(K\)-sweep.

I recommended pushing \(K\) at moderate \(N\) in four consecutive
verdicts (v33, v34, v35, v37). v38 declares the project complete
at \(K = 64\) without running \((1000, 96), (1000, 128), (500, 256)\),
which I predicted would give bounds in the 1.29–1.30 range.

**This is a defensible call**, given that:
* The current bound \(1.2851\) is a clean SOTA exceedance by
  \(0.0051\).
* Pushing \(K\) further would *probably* improve by \(\sim 0.01\)–\(0.02\),
  but is not *certain* to.
* The compute cost of larger \(K\) is \(O(K^{2})\) or \(O(K^{3})\)
  per LP, making the experiments expensive.
* Declaring early gives a clean publishable result.

But it is also worth noting: **the project is unlikely to be at the
v30-framework's true ceiling**. The data \(N\)-saturation curve
and the \(K\)-sensitivity strongly suggest that bounds in the
1.29 range are accessible. If the proposer is willing to spend
~30–60 minutes of additional compute, the bound could plausibly
reach 1.295.

The proposer's call to ship at 1.2851 is justifiable but
conservative.

### Comment 3. The "open work" list is reasonable.

The four open directions listed in §"Open work" are reasonable
extensions:
1. Larger parameter sweeps — straightforward but expensive.
2. Higher-order Lasserre on \((a, b)\) — likely to give a small
   improvement, similar to other constraint additions.
3. Tighter f-side cell bounds via spline interpolation — could
   accelerate convergence in \(N\), but requires non-trivial
   rigorous analysis of spline error bounds.
4. Leverage v9's even-restricted bound 1.69 — interesting if the
   "approximate evenness" can be made rigorous, but no clear
   path is given.

None of these is a low-hanging fruit. The proposer is honest about
the open status.

### Comment 4. Inheritance is clean.

v38 inherits the full v0–v37 chain. No new bugs.

## 5. Feedback to the proposer

If the proposer wants to ship at \(1.2851\), this is a defensible
publishable result. **Recommended**: write a clean technical paper
documenting:

1. The v14 autocorrelation Fejér construction (the key conceptual
   step).
2. The v30 multi-point Fejér mean refinement.
3. The numerical sweep showing convergence to \(1.2851\).
4. The v29 diagnostic showing that v14 alone is loose by a factor
   of \(\sim 8\), motivating v30's pointwise tests.

If the proposer wants to push the bound higher first, run:
* \((1000, 96, 101)\) — predicted \(\sim 1.290\)
* \((1000, 128, 101)\) — predicted \(\sim 1.295\)
* \((500, 256, 101)\) — predicted \(\sim 1.30\)

before writing the summary.

Either path is valid. The bound \(1.2851\) is the project's
clean SOTA-exceeding result.

---

**Final verdict: VALID** (project summary, no new content; the
declared bound \(C_{6.2}\ge 1.2851\) at \((N, K) = (3000, 64)\)
is the v35 result, exceeding the SOTA \(1.28\) by \(0.0051\). The
summary correctly attributes the construction and history. The
project is shippable at this bound, though pushing \(K\) further
could plausibly reach \(\sim 1.295\)).
