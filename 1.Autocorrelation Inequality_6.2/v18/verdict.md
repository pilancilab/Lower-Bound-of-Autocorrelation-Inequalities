# Verdict — Autocorrelation Inequality 6.2, v18

## 1. Summary

v18 adds the **L³ constraint** \(\int F^3 \le \Omega^2\) (via the
chain \(F^3 \le \Omega^2 F\) and \(\int F = 1\)) lifted via Jensen on
\(x^{3}\) to give \(L_F\sum W_m^{3} \le \Omega^{2}\), where \(W_m\)
are F-side cell averages. The lifted constraint is convex via
\(u_m \ge W_m^{3}\) and \(\omega \ge \Omega^{2}\) auxiliary variables
plus the linear \(L_F\sum u_m \le \omega\).

The numerical bound at \((N, K) = (8, 8)\) is \(\Omega^{(18)} = 1.033\),
*identical* to v14 alone. The L³ constraint is mathematically valid
but does not bind on the v14 Fejér optimum.

## 2. Verdict

**VALID** (rigorous, convex, but trivially satisfied).

The L³ constraint is mathematically correct: every admissible \(f\)
gives a feasible point. The convexity is proper (cubic atom on
nonneg variables, square atom on \(\Omega\), linear sum). The
inheritance from v14 is clean.

But the constraint is **structurally vacuous**, and I can show this
by computing the minimum of \(\sum W_m^{3}\) subject to the
constraints.

## 3. Verification of validity, and vacuousness

### The constraint chain. ✓

\(F \ge 0\), \(\int F = 1\), \(F\le \Omega\) pointwise gives
\(F^{3} \le F\cdot \Omega\cdot \Omega = \Omega^{2}F\), and
integrating: \(\int F^{3}\le \Omega^{2}\). ✓

Jensen on \(x\mapsto x^{3}\) (convex on \([0,\infty)\)) applied to
\(F\) with the cell-uniform measure on \(B_m\):
\((1/L_F)\int_{B_m}F^{3}\ge ((1/L_F)\int_{B_m}F)^{3} = W_m^{3}\),
hence \(\int_{B_m}F^{3}\ge L_F W_m^{3}\). Summing over \(m\):
\(\int F^{3}\ge L_F\sum W_m^{3}\). ✓

Combining: \(L_F\sum W_m^{3}\le \Omega^{2}\). ✓

### Convexity. ✓

* \(u_m \ge W_m^{3}\): epigraph of \(x^{3}\) on \(\{x\ge 0\}\),
  which is the convex set \(\{(W, u): W\ge 0, u\ge W^{3}\}\). ✓
* \(\omega\ge \Omega^{2}\): epigraph of square, convex. ✓
* \(L_F\sum u_m\le \omega\): linear. ✓

All conjunctions are convex.

### Validity at the rank-one truth. ✓

For any admissible \(f\), set \(W_m^{*}=\int_{B_m}F/L_F\),
\(u_m^{*} = (W_m^{*})^{3}\), \(\omega^{*}=(\Omega^{*})^{2}\). All
constraints hold at equality. ✓

### Vacuousness — the constraint reduces to \(\Omega\ge 1\).

The LP variables \(W_m\) satisfy \(\sum L_F W_m = \int F = 1\)
(mass), so \(\sum W_m = 1/L_F = 4N\) (with \(L_F = 1/(4N)\)).

By the **power mean inequality** \(M_3\ge M_1\):
\[
\Big(\frac{\sum W_m^{3}}{4N}\Big)^{1/3} \;\ge\; \frac{\sum W_m}{4N} \;=\; 1,
\]
hence \(\sum W_m^{3}\ge 4N\), with equality iff all \(W_m\) are
equal. Multiplying by \(L_F = 1/(4N)\):
\[
L_F\sum W_m^{3} \;\ge\; 1.
\]

So the constraint \(L_F\sum W_m^{3}\le \Omega^{2}\) becomes
\(1\le \Omega^{2}\), i.e., \(\Omega\ge 1\). **This is the trivial
bound** that v0 and every subsequent iteration already had.

The proposer correctly identifies this in the "Numerical observation":
> the L³ constraint **does not bind** at the LP optimum because the
> LP picks \(W\) close to uniform (which minimizes \(\sum W^{3}\)).

The minimizer of \(\sum W^{3}\) subject to \(\sum W = const\) is
the *uniform* distribution (since \(x^{3}\) is convex), and the
uniform \(W\) gives \(\sum W^{3} = 4N\) exactly. So the LP's
relaxation can always pick \(W\) uniform, achieving \(\Omega = 1\)
in the L³ constraint.

## 4. Comments

### Comment 1. The L³ constraint is structurally trivial.

The L³ chain \(\int F^{3}\le \Omega^{2}\int F = \Omega^{2}\) is the
"trivial" \(L^{3}\)-norm bound for nonneg functions: by Hölder,
\(\int F^{3}\le \|F\|_{\infty}^{2}\|F\|_{1} = \Omega^{2}\). Combined
with Jensen on \(W^{3}\), this just gives the trivial
\(\Omega^{2}\ge 1\), i.e., \(\Omega\ge 1\).

For the L³ constraint to give *more* than the trivial bound, you'd
need an *upper* bound on \(\sum W^{3}\) (not a lower bound), which
would force the LP into a non-uniform \(W\). But upper-bounding
\(\sum W^{3}\) requires *additional information* about the
self-convolution structure of \(F\), which isn't in v18.

### Comment 2. To make the L³ constraint useful, you'd need either:

* A tight lower bound on the *convexity defect* \(\sum W^{3} - 4N\),
  which measures how much \(W\) deviates from uniform. This would
  give a *strict* improvement over \(\Omega = 1\) for non-uniform
  \(W\).
* Or a constraint that *forces* non-uniform \(W\) by some other
  mechanism (e.g., a Fourier-side constraint that the uniform
  \(W\) violates).

The v14 autocorrelation Fejér is exactly the "force non-uniform
\(W\)" mechanism via the Fejér sum. So the v14 framework is the
right approach for tightening; v18's L³ adds nothing on top.

### Comment 3. v18 is in the same equivalence class as v15, v16.

v15 added v7/v11 (W) constraints to v14 — didn't help.
v16 added v13 Bochner Toeplitz to v14 — didn't help.
v18 adds the L³ constraint to v14 — doesn't help.

Pattern: the v14 autocorrelation Fejér is *strictly the binding*
constraint at the LP optimum, and adding "necessary but slack"
constraints doesn't tighten the bound. The path forward is to
either increase \((N, K)\) (v17 parameter sweep) or to find a
constraint that the v14 LP optimum *violates* — i.e., a constraint
that the binding \(W\) does not satisfy.

### Comment 4. Don't confuse "valid constraint" with "useful constraint".

v18 is a clean example of a *valid but vacuous* constraint. The
proposer is correctly distinguishing the two: the proof is rigorous
(no algebra errors), but the constraint doesn't tighten the LP. This
is good practice — at least the vacuousness is identified up front
rather than buried.

But v18 should not be counted as "progress" toward the SOTA — it's
*structurally* a no-op, the same as v15 and v16. The next iteration
should focus on the v14 framework's parameter sweep (v17) or the
Fourier-side higher-degree Lasserre lift, not on more L³-style
constraints.

---

**Final verdict: VALID** (rigorous and convex, but trivially
satisfied — same outcome as v15, v16).
