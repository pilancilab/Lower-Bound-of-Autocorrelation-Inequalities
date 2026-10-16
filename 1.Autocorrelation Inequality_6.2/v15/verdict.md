# Verdict — Autocorrelation Inequality 6.2, v15

## 1. Summary

v15 combines all rigorously valid constraints from v14 (autocorrelation
Fejér) and v7/v11 (centered + shifted (W) windowed integral
inequalities). Each constraint is individually valid, and the
combination is valid for the same reason. The numerical bound at
\((N, K) = (8, 8)\) is \(\Omega^{(15)} = 1.033\) — *identical* to v14
alone. The (W) constraints don't bind on top of v14.

## 2. Verdict

**VALID.** No new mathematical content, just a combination of two
previously-validated constraint sets. The validity proof is
inherited from v14 and v7/v11. The relaxation optimum
\(\Omega^{(15)}\le C_{6.2}\) is established by the same rank-one
feasibility argument.

The proposer's empirical observation that the (W) constraints don't
contribute additional tightening at this scale is consistent: the
binding constraint at the LP optimum is the autocorrelation Fejér
(via the chain \(M[k,k]+M[K+k,K+k]\to u_k\to v_k\)), which is *much*
stronger than the (W) constraints' simple convex sums of \(Q[j,k]\).

## 3. Comments

### Comment 1. The (W) constraints are dominated.

The (W) constraints are *necessary conditions* on the second-order
moments \(Q[j,k]\) of \(p\), via the bilinear convolution structure.
The autocorrelation Fejér constraint is a *much sharper* necessary
condition on the *fourth-order* Fourier coefficients
\(|\hat f(k)|^{4}\), via a higher-order SDP lift. So it's not
surprising that the Fejér constraint dominates at the optimum.

This is *informative*: it tells the proposer that further (W) work
is unlikely to help, and the path forward is to strengthen the
Fejér side (more $K$, finer grid, or higher-order Lasserre on the
moment matrix \(M\)).

### Comment 2. The proposer's v16 plan (higher-degree Lasserre on \(M\)) is the right next move.

If the binding constraint is Fejér, and the Fejér constraint is
itself bounded by the v4 SDP lift's tightness (via
\(M[k,k]\ge a_k^{2}\)), then strengthening the SDP lift on \(M\) —
e.g., extending it from degree 1 to degree 2 in \((a, b)\) — should
help. v16 attempts this; I'll verify it next.

---

**Final verdict: VALID** (no improvement, but no error).
