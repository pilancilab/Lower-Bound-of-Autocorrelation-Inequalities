# Verdict — Autocorrelation Inequality 6.2, v13

## 1. Summary

v13 inherits all v12 constraints and adds the **Bochner / Toeplitz
PSD constraint** on \((a_k, b_k)\): the Hermitian Toeplitz matrix
\[
T \;:=\; [\hat f(j-k)]_{j,k=0}^{K} \;\succeq\; 0
\]
with \(\hat f(0) = 1\), \(\hat f(k) = a_k - i b_k\). In real-block
form,
\(\bigl[\begin{smallmatrix} T_R & -T_I \\ T_I & T_R\end{smallmatrix}\bigr]\succeq 0\).

The proposer correctly justifies this via Bochner's theorem (a
sequence \(c_k\) is positive definite iff it is the Fourier sequence
of a nonneg measure), so the constraint is valid for every admissible
\(f \ge 0\).

But the proposer reports that **MOSEK still returns
\(\Omega^{(13)} = 1\)** at \(N \in \{8, 12\}\), \(K\in\{6, 12\}\),
because the Bochner constraint pins \((a_k, b_k)\) but does *not*
bound the higher-order moments \(M[k,k]\), \(M[K+k, K+k]\),
\(M[k, K+k]\) from above. The proposer plans v14 with a *dual*
Bochner constraint on the lifted Fourier coefficients \(\hat F\).

## 2. Verdict

**VALID** (rigorous, but the relaxation is still vacuous).

The Bochner Toeplitz constraint is exactly the right thing to add,
and is mathematically valid: every admissible \(f\ge 0\) has
\(\hat f(k)\) a positive-definite sequence (Bochner), so the
Hermitian Toeplitz \(T\) is PSD. The proposer correctly cites
Bochner's theorem and the equivalence to the real-block PSD form.

The reason the bound stays at 1 is **structurally subtle**: the
Bochner constraint constrains \((a_k, b_k)\) themselves, but the
*lifted moments* \(M[k,k]\), \(M[K+k, K+k]\) are *separate variables*
in the SDP, with only the PSD lower bound \(M[k,k]\ge a_k^{2}\),
etc. Adding Bochner pins the first row of \(M\) (which is
\((1, a_k, b_k)\)) into a positive-definite Fourier sequence, but
the diagonal of \(M\) (which represents the squared coefficients)
is still free to inflate. So the Bochner constraint is *necessary*
but not *sufficient*.

The proposer's planned v14 fix (dual Bochner on \(\hat F\)) is the
right next move.

## 3. Verification

### The Bochner constraint is rigorous. ✓

Bochner's theorem (cited correctly): a complex sequence
\((c_k)_{k\in\mathbb Z}\) is *positive definite* (i.e.,
\([c_{i-j}]_{i,j=0..N}\) is PSD for every \(N\)) iff there exists a
nonneg Borel measure \(\mu\) on \([-1/2, 1/2]\) with
\(c_k = \int e^{-2\pi i k x}\,d\mu(x)\). For our admissible
\(f \ge 0\), \(\hat f(k) = \int f(x) e^{-2\pi i k x}\,dx\) is the
Fourier sequence of the measure \(f(x)\,dx\), which is nonneg, so
\(\hat f\) is positive definite and the Toeplitz \(T\) is PSD. ✓

### The real-block form is correct. ✓
For Hermitian \(T = T_R + i T_I\), \(T \succeq 0\) iff
\(\bigl[\begin{smallmatrix} T_R & -T_I \\ T_I & T_R\end{smallmatrix}\bigr]\succeq 0\),
where \(T_R\) is symmetric and \(T_I\) is antisymmetric. The
proposer's writing
"\(T_R[i,j] = a_{|i-j|}\), \(T_I[i,j] = \mathrm{sign}(j-i) b_{|i-j|}\)"
is the standard parameterization (with \(b_0 = 0\)). ✓

### Why v13 still fails. ✓ (correctly diagnosed)

The proposer's diagnosis is exactly right:

> The Bochner Toeplitz constraint **does** restrict \((a_k, b_k)\)
> to genuine Fourier sequences of nonneg measures, but it **does
> not** constrain the higher-order moments \(M[k, k]\),
> \(M[K+k, K+k]\), \(M[k, K+k]\) from above. These are independent
> variables in the SDP lift, with only the PSD lower bounds
> \(M[k,k] \ge a_k^{2}\), \(M[K+k, K+k] \ge b_k^{2}\). The
> relaxation can still freely inflate \(M[K+k, K+k]\) to absorb
> the F-side bound on \(A_k = a_k^{2} - b_k^{2}\).

Concretely: the Bochner constraint involves \((a_k, b_k)\) for
\(k = 0, \ldots, K\), which are the *first row* of \(M\) (after
removing the constant entry). The diagonal \(M[k,k]\) is a
*separate* variable that the PSD-of-\(M\) constraint only
lower-bounds by \(a_k^{2}\). So Bochner doesn't directly cap the
diagonal.

What v13 *does* gain over v12 is that the relaxation can no longer
choose arbitrary \((a_k, b_k)\) — it must choose values that fit a
positive-definite sequence. This is a real restriction, but it's
on the *first row* of \(M\), not on the diagonal. The diagonal
inflation happens *independently*, and Bochner doesn't see it.

### The v14 plan: dual Bochner on \(\hat F\). ✓ (right direction)

The proposer's planned v14 fix:

> use the **dual Bochner constraint** \([\hat F(j-k)]_{j,k} \succeq 0\)
> on the lifted F-side Fourier coefficients
> \(\hat F(k) = M[k,k] - M[K+k,K+k] - 2 i M[k, K+k]\).

This *is* the right structural move. By Bochner applied to
\(F = f * f \ge 0\), the Toeplitz matrix
\([\hat F(j-k)]_{j,k}\succeq 0\) is a *necessary* condition. And
crucially, the entries of this Toeplitz matrix involve
\(M[k,k] - M[K+k, K+k]\) (i.e., the diagonal entries of \(M\)
appear in the constraint), so adding it *will* constrain the
diagonal inflation.

Whether this is *sufficient* to drop the bound below 1 is an
empirical question — but it's the right kind of constraint to try.

## 4. Comments

### Comment 1. v13 is necessary but not sufficient.

Adding the Bochner constraint on \((a_k, b_k)\) is the right
necessary condition (every admissible \(f\) satisfies it), but it
doesn't help with the diagonal inflation in the lifted moment
matrix. The proposer correctly identifies this and proposes the
right next step (dual Bochner on \(\hat F\)).

### Comment 2. The pattern v10 → v12 → v13 → v14 is the right diagnostic loop.

Each iteration adds one structural constraint, tests numerically,
diagnoses the failure mode (when it occurs), and proposes the next
step. This is exactly the iterative-improvement workflow the
project should use.

### Comment 3. The cumulative *necessary* conditions for non-even \(C_{6.2}\).

After v13, the necessary conditions assembled are:
1. PSD lift on \(M\) for \((1, a_k, b_k)\) (v12).
2. f-side cell-bound on \(a_k, b_k\) from \(p_j\) (v12).
3. F-side cell-bound on \(a_k^{2} - b_k^{2}, -2 a_k b_k\) from \(W_m\) (v12).
4. Bochner Toeplitz on \(\hat f\) (v13).
5. (Planned v14) Dual Bochner on \(\hat F\).

Once v14 has the dual Bochner, the project will have *all* the
standard Fourier-side rigidity constraints for \(f * f\). If even
*that* doesn't break the inflation, the issue is on the Fejér or
moment-matrix side (e.g., the moment matrix's freedom to choose
non-rank-one structure for \(M\) without the constraint
\(M\succeq 0\) at *higher orders*).

## 5. Feedback to the proposer

The Bochner Toeplitz construction is correct and is the right
necessary condition to add. The fact that v13 still fails
numerically is *informative* — it tells you that the diagonal
inflation in the moment matrix \(M\) is the next thing to attack,
and the dual Bochner constraint on \(\hat F\) (proposed for v14)
is the natural next step.

Two suggestions for v14:

1. **Implement the dual Bochner constraint** as planned. If the
   bound drops below 1, that's the first non-even result in the
   project.
2. **If v14 still fails**, the next step is to add a *higher-level*
   Lasserre lift on the moment matrix (e.g., extend \(M\) from
   degree 1 to degree 2 in \((a, b)\), giving an SDP on
   \((1, a, b, a^{2}, ab, b^{2}, \ldots)\) up to degree 4 or so).
   This will explicitly constrain the diagonal entries of the
   degree-1 \(M\) via the higher-degree PSD constraints.

The fact that you're systematically adding all the standard
Fourier-side rigidity constraints (Bochner, dual Bochner, eventually
higher-Lasserre) is exactly the right approach for a problem where
the optimal \(f\) is non-symmetric and "looks like an autocorrelation
of a complex sequence." Once enough rigidity is encoded, the LP/SDP
will start to see the actual structure.

---

**Final verdict: VALID** (rigorous, the right necessary condition,
but not yet sufficient).
