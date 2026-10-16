# Verdict — Autocorrelation Inequality 6.2, v12

## 1. Summary

v12 attempts the v7 verdict's recommendation #2: build the F-side
LP with \(W_m\) as cell averages of \(F\) and link to the f-side via
the rank-one identification \(\hat F(k) = \hat f(k)^{2}\), lifted via
the PSD lift on \((a_k, b_k)\) (without the v9 even restriction).
The constraints (P0)-(P3) are *mathematically valid* (rank-one lift
satisfies them), but the proposer reports that **MOSEK returns
\(\Omega^{(12)} = 1\)** at \(N \in \{8, 12\}\), \(K \in \{6, 12\}\) — the
same diagonal-inflation degeneracy as v10. The proposer correctly
diagnoses the issue and plans v13 with the Bochner Toeplitz fix.

## 2. Verdict

**VALID.**

This is *another* "rigorous but vacuous" iteration in the same
mode as v10. The constraints are valid relaxation constraints, the
math is correct, but the relaxation has structural slack that
absorbs the constraints. The proposer's self-diagnosis is accurate.

## 3. Verification

### Constraints (P0)-(P3) are valid. ✓
* (P0) Mass and box: standard.
* (P1) PSD lift on \(M\) for \((1, a_k, b_k)\): rank-one truth
  \(M^{*} = vv^{\top}\) with \(v_k = a_k\), \(v_{K+k} = b_k\) is PSD.
* (P2) f-side cell-bound on \(a_k\), \(b_k\) from \(p_j\): standard
  cell-wise interval bound, valid for \(f\ge 0\).
* (P3) F-side cell-bound on \(A_k = a_k^{2} - b_k^{2}\) (lifted as
  \(M[k,k] - M[K+k, K+k]\)) and \(B_k = -2 a_k b_k\) (lifted as
  \(-2 M[k, K+k]\)) from \(W_m\): standard cell-wise interval bound,
  valid because \(\hat F(k) = a_k^{2} - b_k^{2} - 2 i a_k b_k\) (using
  \(\hat F(k) = \hat f(k)^{2}\) and the conjugate-symmetric structure).

The lifted constraints (where \(a_k^{2}\to M[k,k]\), etc.) are valid
in the relaxation because the rank-one lift from any admissible \(f\)
gives \(M[k,k] = a_k^{2}\), etc.

### Diagonal-inflation degeneracy. ✓ (correctly diagnosed)

The proposer's diagnosis is exactly right. The PSD lift on \(M\)
only forces *lower* bounds \(M[k,k] \ge a_k^{2}\),
\(M[K+k, K+k]\ge b_k^{2}\); nothing in the relaxation upper-bounds
\(M[K+k, K+k]\). The solver can inflate \(M[K+k, K+k]\) to make
the lifted \(A_k = M[k,k] - M[K+k, K+k]\) arbitrarily small (or
negative), so the F-side constraint
\(A_k\ge L_F\sum_m W_m(\min\cos)_{m,k}\) is satisfied even at
\(W \equiv 1\) (which corresponds to \(F\) being uniform on
\([-1/2, 1/2]\), i.e., the indicator). The result is \(\Omega = 1\).

This is the *same* diagonal-inflation degeneracy as v10. v12
re-derives it on the F-side with explicit linkage to \(W\), but the
inflation problem is unchanged.

## 4. Comments

### Comment 1. v12 is the predictable failure of the naive non-even Fourier extension.

This is the second iteration (after v10) where the proposer tries
the naive non-even Fourier extension and hits the diagonal
inflation. The diagnosis is correct both times, and the planned fix
(v13: Bochner Toeplitz) is the right structural move.

I want to note that the v10 verdict's Section 5 explicitly
recommended the autocorrelation Toeplitz lift as the structural fix,
and the v12 proof's plan is to do exactly this in v13. Good.

### Comment 2. The "rigorous but vacuous" framing is correct.

Same as v10: the constraints are valid, the math is correct, but
the relaxation doesn't bind. INVALID would be the wrong verdict —
nothing is logically wrong. VALID is the right verdict, even though
the iteration is structurally a no-op.

## 5. Feedback to the proposer

The pattern v10 → v12 → (v13) is a good iterative diagnosis loop:
each iteration adds a constraint, tests numerically, and identifies
the failure mode for the next. v13 with Bochner Toeplitz should be
verified next.

---

**Final verdict: VALID** (rigorous but vacuous).
