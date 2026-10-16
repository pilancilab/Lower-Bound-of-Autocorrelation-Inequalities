# Verdict — Autocorrelation Inequality 6.2, v47

## 1. Summary

v47 adds **trigonometric Cauchy–Schwarz upper bounds** on the
diagonal entries of the moment matrix \(M\), coupling them to the
auxiliary \(a_{2k}\) variables introduced in v46:

* **(C1)**: \(M[k, k] \le (1 + a_{2k})/2\) for \(k = 1, \dots, K\).
* **(C2)**: \(M[K+k, K+k] \le (1 - a_{2k})/2\) for \(k = 1, \dots, K\).
* **(C3)**: \(4\,M[k, K+k]^{2}\le 1 - a_{2k}^{2}\) — but this is
  **redundant** given (C1), (C2), and \(M\succeq 0\) (via Schur).

The construction is v47 = v45 + v46 + (C1) + (C2). No numerical
results yet — proposer awaiting theorist verdict.

## 2. Verdict

**VALID** (rigorous, well-motivated tightening of the moment-lift
relaxation).

## 3. Verification

### The Cauchy–Schwarz step. ✓

For \(f\ge 0\) with \(\int f = 1\), define the inner product
\(\langle g, h\rangle_f = \int f\,g\,h\,dx\), which is a positive
semi-definite bilinear form. Cauchy–Schwarz:
\[
\langle g, h\rangle_f^{2} \;\le\; \langle g, g\rangle_f\cdot \langle h, h\rangle_f.
\]
Apply with \(g = 1\), \(h = \cos(2\pi k x)\):
* \(\langle 1, \cos\rangle_f = a_k\)
* \(\langle 1, 1\rangle_f = \int f = 1\)
* \(\langle \cos, \cos\rangle_f = \int f\cos^{2}(2\pi k x)\,dx
  = \int f\cdot \tfrac{1 + \cos(4\pi k x)}{2}\,dx
  = \tfrac{1 + a_{2k}}{2}\) ✓

So \(a_k^{2}\le (1 + a_{2k})/2\). Analogously
\(b_k^{2}\le (1 - a_{2k})/2\). ✓

### Lifted constraint validity.

At any rank-one truth from an admissible \(f\):
* \(M[k, k] = a_k^{2}\) (by definition of the rank-one lift).
* The Cauchy–Schwarz inequality gives \(a_k^{2}\le (1 + a_{2k})/2\),
  hence \(M[k, k]\le (1 + a_{2k})/2\) at the rank-one truth. ✓

Adding (C1) as an LP constraint preserves the rank-one feasible
set, so v47 is a valid relaxation. As a min LP, adding constraints
can only **increase** the optimum:
\[
\Omega^{(47)}\ge \Omega^{(46)} \;(\ge \Omega^{(45)} \;\ge \Omega^{(44)}).
\]

### (C3) is redundant. ✓

The proposer correctly notes:
\[
M[k, K+k]^{2}\le M[k,k]\cdot M[K+k, K+k]\quad\text{(Schur from $M\succeq 0$)}
\]
\[
\le \tfrac{1 + a_{2k}}{2}\cdot \tfrac{1 - a_{2k}}{2} = \tfrac{1 - a_{2k}^{2}}{4}\quad\text{(C1, C2)}.
\]
So (C3) is implied by (C1) + (C2) + \(M\succeq 0\) and need not be
added explicitly. ✓

### The trigonometric identities. ✓

Used:
* \(\cos^{2}\theta = (1 + \cos(2\theta))/2\) ✓
* \(\sin^{2}\theta = (1 - \cos(2\theta))/2\) ✓
* (and \(\cos\theta\sin\theta = \sin(2\theta)/2\) for the analogous
  bound on \(M[k, K+k]\), via \(b_{2k}\) — but this is also
  redundant with C1+C2+Schur).

The substitution converting \(\int f\cos^{2}(2\pi k x)\) into a
linear function of the auxiliary \(a_{2k}\) is standard and
correct.

## 4. Comments

### Comment 1. Why the C–S bounds are non-trivial.

In the SDP without (C1), \(M[k, k]\) is constrained from below by
\(a_k^{2}\) (via Schur on \(\begin{pmatrix} 1 & a_k \\ a_k & M[k,k]\end{pmatrix}\succeq 0\))
but **not** from above. The LP can drive \(M[k, k]\) arbitrarily
large to manipulate \(M[k, k] - M[K+k, K+k]\) (which appears in
the Fejér–Riesz \(\tilde c_k\)) without paying any cost.

(C1) caps \(M[k, k]\) at \((1 + a_{2k})/2 \le 1\), and this cap is
*data-dependent* via \(a_{2k}\). At the rank-one truth, the cap
is tight (\(M[k, k] = a_k^{2} = \) C–S bound, with equality iff
\(\cos\) is a constant multiple of \(1\) on \(\text{supp}(f)\),
which is generically not the case). For typical \(f\), the cap
gives a strict upper bound that the LP must respect.

This is the **first constraint in the project that ties the moment
lift's diagonal to the cell-bound variables \(a_{2k}\)** rather
than to the rank-one quadratic structure alone.

### Comment 2. The CS bounds are sharp.

Equality in C–S holds iff \(\cos\) and \(1\) are linearly
dependent in \(L^{2}(f)\), i.e., \(\cos(2\pi k x)\) is constant
on \(\text{supp}(f)\). This is rare but happens at \(k = 0\) and
edge cases. For generic \(f\), the C–S bound is strict, providing
real tightening.

### Comment 3. Compute scaling.

(C1) and (C2) are simple linear constraints, adding \(2K\)
inequalities to the LP. This is negligible compared to the v45/v46
PSD blocks. So v47 should run at the same speed as v46.

### Comment 4. The auxiliary \(a_{2k}\) variables are essential.

Without v46's extension of the \(a, b\) variables to indices
\(K+1, \dots, 2K\), the constraint (C1) at \(k = K\) would need
\(a_{2K}\), which doesn't exist in v45. v46 + v47 form a natural
pair: v46 introduces the auxiliaries, v47 uses them.

### Comment 5. Inheritance is clean.

v47 inherits from v45 + v46. No new bugs.

## 5. Feedback to the proposer

For v48:

1. **Run v47 at the same parameters as v45/v46's sanity table**:
   confirm \(\Omega^{(47)}\ge \Omega^{(45)}\) at every point.
2. **Run v47 at production parameters**: \((1000, 64)\),
   \((1000, 96)\), \((1000, 128)\), \((500, 256)\). Tabulate
   gain over v45 and v44.
3. **The v44 ceiling at \((1000, 96)\) is 1.2870**. If v47
   pushes this to \(\ge 1.290\), v47 is the new best framework.
4. **Predict**: at \((1000, 96)\), the C–S constraint should
   tighten by \(\sim 10^{-3}\)–\(10^{-2}\) because the LP was
   exploiting the unbounded \(M[k, k]\) freedom.

If v47 + v44 + v45 + v46 collectively reach \(\Omega \ge 1.30\) at
moderate parameters, the project breaks new ground.

---

**Final verdict: VALID** (rigorous trigonometric Cauchy–Schwarz
upper bounds on \(M[k, k], M[K+k, K+k]\) coupling the moment lift
to the auxiliary \(a_{2k}\) variables. \((C_3)\) is redundant
given (C1)+(C2)+Schur. \(\Omega^{(47)}\ge \Omega^{(46)}\) provably;
production sweep pending. This is the first project constraint
that ties \(M[k, k]\) directly to higher-order Fourier coefficients
of \(f\)).
