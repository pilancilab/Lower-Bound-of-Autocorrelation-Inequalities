# Verdict — Autocorrelation Inequality 6.2, v16

## 1. Summary

v16 combines v14 (autocorrelation Fejér) + v13 (Bochner Toeplitz on
\(\hat f\)) + v11 (shifted (W) windows). All constraints are
individually valid. The numerical bound at \((N, K) = (8, 8)\) is
\(\Omega^{(16)} = 1.033\) — *identical* to v14 alone. None of the
added constraints (Bochner or (W)) tightens the v14 Fejér bound at
this scale.

## 2. Verdict

**VALID.** Same situation as v15: the combination of valid
constraints is itself valid, and the LP optimum is a valid lower
bound on \(C_{6.2}\). No new mathematical content; the math is the
union of previously verified material.

The empirical observation is informative: the binding constraint is
the v14 Fejér, and neither v11 (W)s nor v13 Bochner contribute
marginally. This narrows down where to look for the next improvement.

## 3. Comments

### Comment 1. The proposer's diagnosis is correct.

The proof says:
> the binding constraint is the autocorrelation Fejér ...
> the bottleneck is the linear bounds on \(a_k, b_k\) from \(p\)
> (which the LP can satisfy by choosing a non-uniform \(p\) that
> minimizes \(\sum |\hat f(k)|^{4}\)).

This is exactly right. The LP optimum is determined by:
1. The Fejér constraint \(\Omega \ge 1 + 2\sum (1 - k/(K+1))v_k\).
2. The lift \(v_k \ge u_k^{2}\), \(u_k = M[k,k]+M[K+k,K+k] \ge a_k^{2}+b_k^{2}\).
3. The cell-bound \(\sum_j p_j(\min\cos)_{j,k}\le a_k\le \sum_j p_j(\max\cos)_{j,k}\).
4. \(\sum_j p_j = 1\).

The LP minimizes \(\Omega\) subject to (1)–(4), which means it
chooses \(p\) (and via (3), the \(a_k, b_k\)) to minimize the Fejér
sum. The Bochner Toeplitz on \((a_k, b_k)\) and the (W) constraints
on \(Q[j,k]\) are *necessary* conditions on \(f\), but they don't
directly enter the Fejér chain — they only restrict \(p, Q\), which
the LP can already adjust freely within those constraints.

### Comment 2. The path forward is to strengthen the \(p\)-to-\(a_k\) link.

If the bottleneck is the cell-bound on \(a_k\) from \(p\), then the
right move is to make this link tighter. Two options:

1. **Higher-order moment matching.** Add constraints
   \(\sum_j p_j x_j^{m}(\min\cos)\le \int x^{m} f \cos(2\pi k x)\le \sum p_j x_j^{m}(\max\cos)\)
   for \(m = 1, 2, \ldots\), pinning down higher-order moments of
   \(f\). This forces \(p\) to be a more accurate discretization of
   the true \(f\).
2. **Quadrature error tightening.** The current cell bounds use the
   worst-case min/max of \(\cos(2\pi k x)\) on \(I_j\). For large
   \(k\), this is loose because the cosine oscillates within a cell.
   Use a midpoint or trapezoidal rule with explicit error bound, as
   in 6.3 v1's Lemma 6.

Either of these would tighten the \(p\to (a_k, b_k)\) link and let
the Fejér constraint actually bind on a more constrained set.

### Comment 3. v17 (parameter sweep) is also a reasonable next step.
The proposer's plan to run v17 at large \((N, K)\) is informative: if
the bound saturates at some value below the asymptotic \(4/3\), that
tells us the discretization is the bottleneck.

---

**Final verdict: VALID** (no improvement, no error).
