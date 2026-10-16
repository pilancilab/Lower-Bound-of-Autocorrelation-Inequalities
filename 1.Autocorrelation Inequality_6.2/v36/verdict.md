# Verdict — Autocorrelation Inequality 6.2, v36

## 1. Summary

v36 takes v30 (autocorrelation Fejér + pointwise Fejér mean tests
on \(F\)) and adds the v28 Bochner Toeplitz PSD constraint on
\(F\)'s Fourier sequence (the Hermitian Toeplitz of the lifted
\(\hat F(k) = (a_k - i b_k)^{2}\) values).

Numerical result at \((N, K, n_{\text{tests}}) = (16, 16, 33)\):
\(\Omega^{(36)} = 1.147\), **identical** to v30 at the same
parameters. The Bochner constraint is not binding.

## 2. Verdict

**VALID** (rigorous, but the new constraint produces no
improvement — same observation as v28 when added on top of v14).

## 3. Verification

Both ingredients were previously verified:
* v30 pointwise Fejér mean tests + v14 SOC chain — verified in
  `v30/verdict.md`.
* v28 Bochner Toeplitz PSD on \(\hat F\) — verified in
  `v28/verdict.md` (Bochner's theorem applied to \(F = f * f\ge 0\)).

The conjunction is valid: the rank-one truth from any admissible
\(f\) satisfies all constraints simultaneously, so the LP feasible
set still contains the rank-one truth and \(\Omega^{(36)}\le C_{6.2}\). ✓

The numerical observation \(\Omega^{(36)} = \Omega^{(30)} = 1.147\)
at \((16, 16)\) is consistent with the v28 finding: the SOC chain
\(v_k\ge u_k^{2}\) already forces the LP to behave as if rank-one
at the optimum, and rank-one truths automatically satisfy Bochner.
So adding Bochner on top of v30 excludes nothing the LP would
have picked anyway.

## 4. Comments

### Comment 1. Same observation as v28, on a different base.

The v28 verdict already established that adding Bochner Toeplitz
on \(\hat F\) on top of v14 is vacuous, for the structural reason
that v14's SOC chain forces rank-one behavior. v36 reproduces this
observation on top of v30 (which is v14 + extra constraints, all of
which preserve rank-one).

This is **expected** and adds no new information about the
framework's reach.

### Comment 2. Best bound is unchanged.

The best rigorous bound on \(C_{6.2}\) remains v35's
\(1.2851\) at \((N, K) = (3000, 64)\). v36 does not change it.

The proposer correctly cites the v35 result in the boxed claim.

### Comment 3. Inheritance is clean.

v36 inherits from v30 + v28. No new bugs.

### Comment 4. The proposer should stop adding orthogonal constraints to v30.

The v15–v28, v32, v36 history all confirm that adding further
constraint families on top of the core v14/v30 framework is
**always vacuous**. The reason: the SOC chain forces the LP
optimum to lie at a rank-one face, and rank-one truths
automatically satisfy any nonneg-functional constraint (Bochner,
Toeplitz, moment-cone, etc.).

The only ways forward are:
1. **Larger parameters** \((N, K)\) — what v33/v34/v35 are doing.
2. **A fundamentally different relaxation** that breaks the
   rank-one boundary — e.g., higher-order Lasserre on \(f\) (not
   just on \((a, b, M)\)), or a dual SDP approach, or a SOS
   hierarchy on the variational problem of \(f\).

Adding more constraints to v30 will not help; it has been verified
empirically across nearly 20 iterations.

## 5. Feedback to the proposer

For v37: **stop adding constraints to v30** and run the
\(K\)-sweep at moderate \(N\):

* \((1000, 96)\), \((1000, 128)\), \((2000, 96)\), \((2000, 128)\),
  \((500, 256)\)

Each of these is more likely to give a real improvement than any
new constraint type. The expected bounds are in the 1.290–1.300
range.

If the bound saturates somewhere in the 1.29–1.30 region, document
the asymptotic curve and declare the v30 framework's ceiling.

---

**Final verdict: VALID** (rigorous, but vacuous; the Bochner
Toeplitz constraint on \(\hat F\) is not binding when added on
top of v30, just as it was not binding when added on top of v14
in iteration v28. Best bound on \(C_{6.2}\) remains v35's
\(1.2851\)).
