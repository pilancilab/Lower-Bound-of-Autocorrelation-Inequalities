# Verdict — Autocorrelation Inequality 6.2, v46

## 1. Summary

v46 extends v45's Bochner Toeplitz constraint from degree \(K\) to
degree \(2K\). The proposer adds auxiliary LP variables
\(a_{K+1}, \dots, a_{2K}\) and \(b_{K+1}, \dots, b_{2K}\), with
their own cell-wise rectangular bounds from \(p\), and enforces the
larger \((2K+1)\times(2K+1)\) Hermitian Toeplitz PSD constraint
\(T_{2K} \succeq 0\) on the full sequence
\((c_0, c_1, \dots, c_{2K})\) with \(c_k = a_k - i b_k\).

The Fejér–Riesz constraint and the v14 SOC chain still act on
\(a_{1..K}, b_{1..K}\) only.

**No numerical results yet** — proposer awaiting theorist verdict.

## 2. Verdict

**VALID** as a relaxation; mathematically sound.

## 3. Verification

### (1) Cell bounds on the auxiliary variables. ✓

For each \(k\in\{K+1, \dots, 2K\}\), the auxiliary
\(a_k = \int f \cos(2\pi k x)\,dx\) (and \(b_k\) analogously) is a
linear functional of \(f\). Per cell \(j\), the maximum (resp. minimum)
of \(\cos(2\pi k x)\) over the cell gives a valid pointwise bound on
the per-cell contribution \(p_j\cdot \cos\), so the cell-wise
rectangular bound
\[
\sum_j p_j \min_{x\in C_j}\cos(2\pi k x) \;\le\; a_k \;\le\; \sum_j p_j \max_{x\in C_j}\cos(2\pi k x)
\]
is valid for any nonneg \(f\) with cell mass \(p_j\). Same for \(b_k\). ✓

### (2) Bochner–Toeplitz at degree \(2K\). ✓

By Bochner's theorem applied to the positive measure \(d\mu = f\,dx\)
on the torus, the Hermitian Toeplitz matrix
\(T_{2K}[i, j] = c_{i-j}\) for \(0\le i, j\le 2K\) is positive
semidefinite. The proof is identical to v45's verification
(\(z^{*}T_{2K}z = \int f|\phi|^{2}\,dx \ge 0\)), just at a larger
matrix size.

The real-lift encoding \(\widetilde T_{2K} = \begin{pmatrix} T_R & -T_I \\ T_I & T_R\end{pmatrix}\succeq 0\)
is the same standard equivalence as v45. ✓

### (3) Validity of the relaxation.

For every admissible \(f\), the rank-one truth in v45 extends
naturally to v46: the actual Fourier coefficients
\(a_k, b_k\) for \(k = K+1, \dots, 2K\) satisfy the cell bounds
(with equality at the rank-one truth's specific cell densities),
and the full \(T_{2K}\) is PSD by Bochner. So every rank-one truth
remains feasible, and v46 is a valid relaxation. ✓

### (4) Strict tightness over v45.

v45's constraint is "the truncated Toeplitz \(T_K\) (\((K+1)\times(K+1)\)) is PSD."
v46's constraint is "there exist \(a_{K+1..2K}, b_{K+1..2K}\) in
their cell-bound boxes such that the full \(T_{2K}\) is PSD." This
is strictly stronger, because:

* If \((a_1, \dots, a_K, b_1, \dots, b_K)\) admits a PSD extension
  via some \(a_{K+1..2K}\), it also satisfies the truncated PSD
  constraint at degree \(K\) (the truncated principal submatrix of
  a PSD matrix is PSD).
* The converse fails: a PSD truncated Toeplitz \(T_K\) need not
  admit any PSD extension to \(T_{2K}\), let alone one with the
  cell bounds. (Standard moment-problem fact: not every truncated
  PSD moment sequence extends to an infinite positive measure.)

So **v46's constraint is strictly stronger than v45's**, and the
LP optimum can only increase (or hold).

\(\Omega^{(46)} \ge \Omega^{(45)}\). ✓

### (5) The auxiliary variables don't relax the LP.

Adding new LP variables would normally **relax** a minimization LP
(more freedom for the optimizer). But here the new variables are
**only** referenced in the new constraint (Bochner at degree \(2K\)),
not in any constraint that involves the original variables. So the
new variables can be freely chosen by the LP to satisfy the
Bochner constraint, and:

* If for some \((a_1, \dots, a_K)\) the larger Bochner is satisfiable
  with some choice of \(a_{K+1..2K}\), the LP picks that and the
  combined point is feasible.
* If for some \((a_1, \dots, a_K)\) no such choice exists, the LP
  cannot pick \((a_1, \dots, a_K)\) at all.

So the **net effect** of v46 over v45 is to constrain the original
\((a_1, \dots, a_K)\) more tightly (must admit a PSD extension), not
to add freedom. ✓

## 4. Comments

### Comment 1. v46 follows the standard moment-extension pattern.

The truncated moment problem (Hausdorff/Carleman/Hamburger) is
classical: for a sequence to be the moments of a positive measure,
all truncated moment matrices must be PSD. **Not all truncated
PSD sequences extend** — there are sequences that satisfy the
truncated PSD condition at degree \(K\) but not at degree \(K + 1\).
v46's degree-\(2K\) constraint exploits this gap to tighten v45.

This is a standard SOS/moment-hierarchy refinement (Lasserre level
\(K\) → level \(2K\) on the trigonometric moment problem).

### Comment 2. The auxiliary variables' cell bounds are loose.

The cell-wise rectangular bounds \([\min_{x\in C_j}\cos(2\pi k x), \max_{x\in C_j}\cos(2\pi k x)]\)
become very loose at large \(k\), because \(\cos(2\pi k x)\) oscillates
rapidly within each cell. At \(k = 2K\) and cell width \(L_f = 1/(4N)\),
the cell can contain \(\sim L_f \cdot 2K = K/(2N)\) full periods of
\(\cos\), so the min and max are \(\pm 1\) and the bound is just
\(\sum_j p_j\cdot 1 = 1\), \(\sum_j p_j\cdot(-1) = -1\) — the trivial
bound \(|a_k|\le 1\).

So at large \(k/N\), the cell bounds are essentially trivial. The
LP relies on the Bochner extension constraint (which couples the
new variables to the old ones) to do the tightening, not on the
cell bounds.

This is fine for validity but means the real tightening comes
from Bochner's coupling effect, not from the cell bounds.

### Comment 3. Compute scaling.

The new \((2K+1)\times(2K+1)\) Bochner block is \(\sim 4\times\)
the v45 block. At \(K = 96\), that's \(\sim 4\times 9000 = 36000\)
LMI entries. MOSEK should handle this in 1–2× the v45 solve time.

### Comment 4. Expected gains.

The proposer wisely says "hard to predict." My estimate:
* At small \(K\) (\(K = 32\)), the truncated moment problem is
  often nearly self-extending, so the gain is \(\sim 10^{-5}\).
* At larger \(K\) (\(K = 64, 96\)), the gap between truncated PSD
  and extendable PSD grows, and the gain may be \(\sim 10^{-4}\)–\(10^{-3}\).
* The combined v44 + v45 + v46 sequence should give
  \(\sim 10^{-3}\) total improvement over v44 at \(K = 96\).

But this is speculative; the proposer should run and report.

### Comment 5. Inheritance is clean.

v46 inherits from v45. No new bugs.

## 5. Feedback to the proposer

**Math is sound, please run.** For v47:

1. **Run v46 at the same parameters as v45's sanity table**:
   \((8, 8), (16, 16), (16, 32), (24, 24), (32, 32), (50, 32)\),
   and verify \(\Omega^{(46)} \ge \Omega^{(45)}\).
2. **Run v46 at production parameters**: \((1000, 64)\),
   \((1000, 96)\), \((1000, 128)\), \((500, 256)\). Tabulate
   \(\Omega^{(46)} - \Omega^{(45)}\) at each.
3. **Predict whether v46 reaches 1.290+** at any \((N, K)\).
4. **If v46 saturates at the same level as v45**, the next step
   is v47 = degree-\(3K\) Bochner (third level of the hierarchy).

The current best bound on \(C_{6.2}\) remains v44's \(1.2870\)
at \((1000, 96)\); v45 and v46 should match or exceed this once
the production sweeps complete.

---

**Final verdict: VALID** (rigorous extension of v45's Bochner
constraint from degree \(K\) to degree \(2K\); the cell bounds, the
Bochner Toeplitz at the larger size, and the standard moment-extension
argument are all correct. \(\Omega^{(46)}\ge \Omega^{(45)}\)
provably; production sweep pending).
