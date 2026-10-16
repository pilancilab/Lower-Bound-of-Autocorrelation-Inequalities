# Verdict — Autocorrelation Inequality 6.2, v11

## 1. Summary

v11 inherits all of v9 (Fejér-kernel test for the even-restricted
case, VALID per v9/verdict.md) and adds **shifted (W) windows**
\(B = [c-h, c+h]\) at offsets \(c\in\{0,\pm L,\ldots,\pm 1/4\}\) and
the same set of half-widths \(h\in\{2L, 3L, \ldots, 1/2\}\) used in
v7. The numerical bound at \((N, K) = (8, 16)\) improves from
v9's 1.349 to v11's **1.376** — a modest but real tightening.

## 2. Verdict

**VALID.**

The v7 (W) claim
\[
\sum_{(j,k)\in\mathcal C(B)}Q[j,k]\le h_B\,\Omega
\]
holds for *any* window \(B = [a,b]\subset[-1/2,1/2]\), not just
centered. The proof in v7 uses Fubini on cell rectangles + the
"fully contained cells lower-bound the iterated integral" + the
trivial \(\int_B F\le h_B\Omega\), none of which require \(B\) to be
centered at zero. So shifted windows give valid constraints for
the same reason centered ones do, and the rank-one lift from any
admissible even \(f\) satisfies them.

The cell-membership criterion \((c-h+1/2)/L + 2\le j+k\le (c+h+1/2)/L\)
is the same as v7's, just instantiated at general \((a, b) = (c-h, c+h)\).

## 3. Numerical observation

v11 reports \(\Omega^{(11)} = 1.376\) at \((N, K) = (8, 16)\),
compared to v9's 1.349 at the same parameters. The improvement is
\(\approx 0.027\), which is real but modest. The proposer notes that
larger \((N, K)\) runs are computationally intensive and predicts a
larger improvement at finer resolutions, since the additional
windows resolve fine-scale features that the centered ones don't see.

For comparison, v9's best is **1.689 at \((50, 64)\)**, which
remains the best rigorous bound on \(C^{\text{even}}_{6.2}\) in the
project. v11 hasn't yet been run at the same large \((N, K)\), so
the actual improvement at scale is unknown.

## 4. Comments

### Comment 1. v11 is the right cumulative move on the moment side.

Both centered and shifted windows are valid (W) constraints, and
adding them all (rather than just centered) gives a strictly tighter
relaxation. This is exactly the v7 verdict's Section 5
recommendation #1 ("add more windows").

### Comment 2. Still bounds \(C^{\text{even}}_{6.2}\), not \(C_{6.2}\).

Like v7 and v9, v11 only handles the even-restricted case. The v12
plan (next iteration) is to attempt the F-side autocorrelation
Toeplitz lift recommended in v7 verdict feedback item 2. I'll
verify v12 separately.

### Comment 3. Inheritance is clean.
v11 inherits from v9 (which inherits from v7), and v7 dropped the
v0 hand-wave. So v11 has no inherited rigor issues.

## 5. Feedback to the proposer

This is a clean cumulative tightening of v9. Two suggestions:

1. **Run v11 at \(N = 50, K = 64\)** (the same parameters as v9's
   best) to see how much the shifted windows actually buy at
   scale. If the improvement is significant (say, above 1.7),
   that's a real numerical advance.
2. **The non-even extension** (v12+) is the bigger structural
   move. Don't get stuck adding more windows to the even case;
   v9 + v11 shifted windows is probably about as much as you can
   squeeze out of the (W) family, and the marginal returns will
   diminish.

---

**Final verdict: VALID.**
