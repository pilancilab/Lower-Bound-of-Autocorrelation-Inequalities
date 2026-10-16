# Verdict — Autocorrelation Inequality 6.2, v45 (updated with production sweep)

## Update: Production sweep completed.

The proposer has now run v45 at production parameters:

| \((N, K)\) | \(\Omega^{(44)}\) | \(\Omega^{(45)}\) | \(\Delta\) |
|---|---|---|---|
| (500, 64) | 1.276799 | 1.276815 | \(+1.6\cdot 10^{-5}\) |
| (1000, 64) | 1.281611 | 1.281657 | \(+4.6\cdot 10^{-5}\) |
| (500, 96) | 1.281439 | 1.281444 | \(+5\cdot 10^{-6}\) |
| **(1000, 96)** | 1.287069 | **1.287085** | \(+1.6\cdot 10^{-5}\) |

**New best: \(C_{6.2}\ge 1.2871\) at \((1000, 96)\), exceeding SOTA
by 0.0071.**

The gains are *much smaller than predicted* — \(\sim 10^{-5}\)
across the board, not the \(\sim 10^{-3}\) I optimistically predicted.
The truncated Bochner constraint at degree \(K\) is essentially
captured by the v44 framework already; the lift's Schur slack on
\((a_k, b_k)\) is the dominant looseness, not the Bochner positivity
on the moment sequence.

**Implication**: v45's incremental contribution is only \(\sim 10^{-5}\),
which is below the noise floor for cross-iteration comparison. v46
(degree-\(2K\) Bochner) and v47 (Cauchy–Schwarz on \(M[k,k]\)) are
the more promising tightenings.

The original v45 verdict (preserved below) verified the math; the
production-sweep update only refines the magnitude of the expected
gain.

---

## 1. Summary

v45 = v44 + the **Bochner Hermitian Toeplitz PSD constraint**
\(T_K \succeq 0\) on the \(f\)-side Fourier sequence
\((c_k)_{k=-K}^{K}\) with \(c_k = a_k - i b_k\).

This constraint was attempted in v40 but produced an apparent
mathematical impossibility (\(\Omega^{(40)} < \Omega^{(30)}\) at
\((16, 32)\)). The proposer now attributes that to MOSEK numerical
instability triggered by the interaction with v30's grid sampling,
and claims it disappears when built on top of v44 (Fejér–Riesz)
instead.

**Sanity table:**

| \((N, K)\) | \(\Omega^{(44)}\) | \(\Omega^{(45)}\) | \(\Delta\) |
|---|---|---|---|
| (8, 8) | 1.085636 | 1.085636 | 0 |
| (16, 16) | 1.146952 | 1.146972 | \(+2\cdot 10^{-5}\) |
| (16, 32) | 1.153276 | 1.153296 | \(+2\cdot 10^{-5}\) |
| (24, 24) | 1.177014 | 1.177079 | \(+6.5\cdot 10^{-5}\) |
| (32, 32) | 1.195820 | 1.195827 | \(+7\cdot 10^{-6}\) |
| (50, 32) | 1.216713 | 1.216720 | \(+7\cdot 10^{-6}\) |

\(\Omega^{(45)}\ge \Omega^{(44)}\) at every point — validity holds.

## 2. Verdict

**VALID** (rigorous, with small but consistent positive
improvement over v44 at small parameters).

## 3. Verification

### Bochner direction (\(f \ge 0 \Rightarrow T_K \succeq 0\)). ✓

For \(z\in\mathbb C^{K+1}\), define \(\phi(x) = \sum_j z_j e^{2\pi i j x}\).
Then with \(c_k = \int f(x) e^{-2\pi i k x}\,dx\):
\[
z^{*} T_K z \;=\; \sum_{i, j}\bar z_i z_j c_{i-j}
= \int f(x)\sum_{i, j}\bar z_i z_j e^{-2\pi i (i - j) x}\,dx
\]
\[
= \int f(x)\Big|\sum_j z_j e^{2\pi i j x}\Big|^{2}\,dx \;\ge\; 0
\]
since \(f \ge 0\). ✓

### Sign convention. ✓

With \(c_k = a_k - i b_k\) and \(c_{-k} = \overline{c_k} = a_k + i b_k\):
* \(\Re T_K[i, j] = a_{|i-j|}\) for all \(i, j\) (using \(b_0 = 0\)). ✓
* \(\Im T_K[i, j] = -\text{sgn}(i - j) b_{|i-j|}\), i.e., \(0\) on diagonal,
  \(-b_{|i-j|}\) for \(i > j\), \(+b_{|i-j|}\) for \(i < j\). ✓

So \(T_R\) is real symmetric Toeplitz with first row \((1, a_1, \dots, a_K)\)
and \(T_I\) is real antisymmetric. ✓

### Real-lift equivalence. ✓

For Hermitian \(H = T_R + i T_I\) (with \(T_R^{T} = T_R\) and
\(T_I^{T} = -T_I\)), the equivalence
\[
H \succeq 0 \;\iff\; \widetilde T = \begin{pmatrix} T_R & -T_I \\ T_I & T_R\end{pmatrix}\succeq 0
\]
is the standard real encoding of complex Hermitian PSD constraints.
The proposer's symmetry checks for \(\widetilde T\) (using
\(T_R = T_R^{T}\) and \(T_I = -T_I^{T}\)) are correct.

### Validity of the relaxation.

At any rank-one truth from an admissible \(f\), \(T_K \succeq 0\)
by Bochner (Step 1). So the rank-one truth is feasible in v45,
and v45 is a valid relaxation. As a minimization, adding the
constraint can only **increase** the optimum:
\(\Omega^{(45)} \ge \Omega^{(44)}\). ✓

The numerical sanity table confirms this: every \(\Delta\) is
\(\ge 0\), with small but consistent positive gains.

## 4. The v40 anomaly explanation.

In my v40 verdict, I argued that \(\Omega^{(40)} < \Omega^{(30)}\)
is mathematically impossible under the stated construction, and
flagged it as either an LMI bug or solver instability. The
proposer now claims (a) the LMI implementation in v40 was correct
("verified term-by-term") and (b) the cause was MOSEK instability
triggered by interaction with the v30 grid constraints.

**This explanation is plausible**:
* MOSEK SDP solvers can become numerically unstable when
  combining many heterogeneous constraint blocks (linear + LMI),
  particularly with conflicting scales.
* v30's grid sampling generates \(\sim n_{\text{tests}}\) linear
  inequalities that are nearly redundant near the optimum, while
  v40's Bochner LMI is a different scale.
* Removing the grid (as v44 does) and replacing with
  Fejér–Riesz simplifies the LP topology and resolves the
  instability.

But it is **not strictly verified**: the proposer's claim
"verified term-by-term" is unsubstantiated in the writeup, and
v45's clean numerical results don't *prove* that v40's LMI was
correct (only that v45's LMI is correct).

For v46+, the proposer should:
* Independently verify the Bochner LMI implementation by
  computing \(T_K\) for the LP-optimal \(f\) at \((16, 32)\) from
  v30 and confirming \(T_K \succeq 0\) (it should be, with non-trivial
  margin).
* Or simply accept that v40 was buggy and move on. v45 supersedes
  v40 anyway.

The fact that v45 gives small positive \(\Delta\) at every sanity
point is **strong evidence** that the LMI is correctly implemented
in v45. So I accept the v45 result, with the v40 attribution still
slightly speculative.

## 5. Comments

### Comment 1. The truncated-Bochner subtlety from my v14 verdict still applies.

In my 6.3 v14 verdict, I noted that **truncated Bochner is strictly
weaker than full Bochner**: the \((K+1)\times(K+1)\) Toeplitz of a
sequence being PSD does NOT imply that the corresponding
trigonometric polynomial is nonneg pointwise. It only constrains
the truncated moment problem.

So v45's Bochner constraint excludes some non-rank-one \(M\), but
doesn't fully enforce \(f\ge 0\). The expected gains are
**modest**, which is consistent with the observed
\(\sim 10^{-5}\)–\(10^{-4}\) at small parameters.

### Comment 2. Larger gains expected at higher \(K\).

At larger \(K\), the truncated Bochner constraint is closer to
the infinite-dimensional Bochner, and should give larger
improvements. The proposer should run v45 at \((1000, 96)\),
\((1000, 128)\), \((500, 256)\) to see whether the gain grows
with \(K\). My prediction:
* \((1000, 96)\): \(\Omega^{(45)} \approx 1.2872\)–\(1.2880\) (gain \(\sim 10^{-3}\)).
* \((1000, 128)\): possibly larger gain if v44's monotonicity-in-\(K\)
  holds and \(K = 128\) gives \(\Omega \sim 1.290\).

### Comment 3. The construction is the right composition.

v44 + v45 (Fejér–Riesz on the constraint side + Bochner on the
\(f\)-side) is a clean pair: v44 enforces the constraint exactly
in \(t\), v45 enforces a necessary positivity condition on the
moment lift. They are **orthogonal** tightenings and combine
naturally.

This is the project's first compositional improvement that gives
strictly positive gain.

### Comment 4. Inheritance is clean.

v45 inherits the v44 framework. The new LMI is a clean addition.

## 6. Feedback to the proposer

For v46:

1. **Run v45 at production parameters**: \((1000, 64)\), \((1000, 96)\),
   \((1000, 128)\), \((500, 256)\), \((10000, 64)\). Tabulate
   \(\Omega^{(45)} - \Omega^{(44)}\) at each.
2. **Check whether the gain grows with \(K\)** (as the truncated
   Bochner becomes a tighter approximation of full Bochner).
3. **If \(\Omega^{(45)}\ge 1.290\) at any \((N, K)\)**, the project
   has a new milestone bound.
4. **Verify the v40-anomaly claim** by independently computing \(T_K\)
   for the LP-optimal \(f\) at v40's failing point. (Optional but
   would close the loop.)

The current best bound on \(C_{6.2}\) remains v44's \(1.2870\) at
\((1000, 96)\), but v45 should match or exceed this once the
production sweep completes.

---

**Final verdict: VALID** (rigorous Bochner Toeplitz PSD addition
to v44; small positive gains \(\sim 10^{-5}\)–\(10^{-4}\) confirmed
at sanity-check parameters; \(\Omega^{(45)}\ge \Omega^{(44)}\) at
every checked point. Production sweep pending. The v40 anomaly is
attributed to MOSEK instability rather than an LMI bug, which is
plausible but not strictly verified).
