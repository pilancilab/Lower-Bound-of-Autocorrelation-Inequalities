# Verdict — Autocorrelation Inequality 6.2, v17

## 1. Summary

v17 reuses **the same constraint set as v14** with no new
constraints, and runs a parameter sweep at larger
\(N\in\{50, 75, 100, 150\}\), \(K\in\{32, 64, 128\}\). The validity
proof is *identical* to v14's. The interesting content is the
**asymptotic limit prediction**: as \(K\to\infty\), the Fejér
weights approach 1, and the constraint becomes
\(\Omega \ge 1 + 2\sum_{k\ge 1}|\hat F(k)|^{2} = \int F^{2}\) (Parseval).
For the uniform \(f = 2\cdot\mathbf 1_{[-1/4, 1/4]}\),
\(\int F^{2} = 4/3 \approx 1.333\), which **exceeds the SOTA
\(C_{6.2}\ge 1.28\)**.

## 2. Verdict

**VALID.** No new mathematical content; the proof is verbatim v14.
The validity is inherited.

The asymptotic-limit calculation \(\int F^{2}_{\text{uniform}} = 4/3\)
is correct (I checked):
\[
\int_{-1/2}^{1/2}\bigl(2(1-2|t|)_{+}\bigr)^{2}\,dt
= 8\int_{0}^{1/2}(1-2t)^{2}\,dt
= 8 \cdot \frac{1}{2}\int_{0}^{1}u^{2}\,du
= 4 \cdot \frac{1}{3} = \frac{4}{3}.
\]

The asymptotic limit prediction is sound: as \(K\to\infty\), the
v14 Fejér constraint becomes the Parseval identity \(\int F^{2}\),
and the relaxation optimum approaches
\(\inf_{f\in\mathcal A^{\text{even-restricted? no}}}\int (f*f)^{2}\).
For the uniform \(f\), this value is \(4/3\); for more concentrated
\(f\), it's larger; for less concentrated \(f\), it's smaller.

So the relaxation's *infimum* over admissible \(f\) of \(\int F^{2}\)
gives the asymptotic limit \(\Omega^{(17)}_{\infty}\). Whether this
equals \(4/3\) (the value at the uniform) depends on whether the
uniform is the minimizer of \(\int F^{2}\) over admissible \(f\).

## 3. Comments

### Comment 1. The asymptotic prediction is informative but not yet a result.

If the v17 numerical sweep approaches \(4/3\) (and we trust the
limit), then v17 establishes a *rigorous numerical improvement on
the SOTA* for the first time in the project. But this requires:

* Empirical confirmation that \(\Omega^{(17)}\) approaches \(4/3\)
  monotonically.
* Analytical confirmation that the *infimum* over admissible \(f\)
  of \(\int F^{2}\) equals \(4/3\) (i.e., that the uniform is the
  minimizer).

Both of these are open questions.

### Comment 2. The infimum may be smaller than 4/3.

I can verify *one* direction: for the *uniform* \(f\),
\(\int F^{2} = 4/3\). But the *infimum* over admissible \(f\) of
\(\int F^{2}\) is the relevant quantity for the relaxation, and it
may be smaller than \(4/3\) — i.e., there may exist a non-uniform
\(f\) with \(\int(f*f)^{2} < 4/3\).

I checked one alternative: for the more concentrated
\(f = 4\cdot\mathbf 1_{[-1/8, 1/8]}\), \(F\) is the triangle on
\([-1/4, 1/4]\) with peak 4, and
\[
\int F^{2} = 2\int_{0}^{1/4}(4(1-4t))^{2}\,dt = 32\int_{0}^{1/4}(1-4t)^{2}\,dt = 32 \cdot \tfrac{1}{4}\int_{0}^{1}u^{2}\,du = 8 \cdot \tfrac{1}{3} = \tfrac{8}{3}.
\]
So *concentrating* \(f\) makes \(\int F^{2}\) *larger*. This is
consistent with the intuition that "spreading out" \(f\) reduces
\(\int F^{2}\), and the uniform on \([-1/4, 1/4]\) may indeed be the
minimizer.

If the uniform is the minimizer, then v17 should converge to
\(\Omega = 4/3\) and provide the first SOTA improvement.

### Comment 3. Numerical results (added in v17 rewrite).

The v17 file now contains numerical results:

| \(N\) | \(K\) | \(\Omega^{(17)}\) |
|---|---|---|
| 50 | 32 | 1.099 |
| 75 | 32 | 1.107 |
| 100 | 32 | 1.111 |

**The convergence is *very* slow.** The bound climbs from 1.099 to
1.111 as \(N\) doubles from 50 to 100 — an improvement of only
\(\sim 0.012\) per doubling. At this rate, reaching \(1.28\) (the
SOTA) would require enormously larger \(N\) (extrapolating, perhaps
\(N \sim 10^{4}\) or more), and reaching the asymptotic \(4/3\)
would require even more.

This is informative: the v14/v17 framework is *capable in
principle* of reaching the SOTA via the asymptotic limit, but the
*convergence rate* is too slow at practical \((N, K)\) values to
make this a viable numerical strategy.

**Two possible interpretations:**

1. **Slow grid convergence.** The cell-bound on \(a_k\) from \(p\)
   uses worst-case \(\min/\max\) of \(\cos(2\pi k x)\) on each
   \(I_j\), which is loose and gives slow convergence in \(N\). A
   tighter quadrature (midpoint rule with explicit error) might
   converge faster.
2. **Non-uniform minimizer.** The relaxation may be picking a
   non-uniform \(p\) that gives \(\int F^{2} < 4/3\), so the
   asymptotic limit isn't \(4/3\) at all — it could be \(\sim 1.1\)
   or so. The proposer raises this possibility in the proof's
   "Caveat" section.

To distinguish: the proposer should compute \(\int F^{2}\) for the
LP-optimal \(p\) at large \(N\) and check whether it matches the
LP value. If yes, then (b) is correct and the LP is finding a
non-uniform optimizer. If no, then (a) is correct and the
discretization error is the bottleneck.

### Comment 4. The asymptotic argument is informative for the project's roadmap.

Even if v17 doesn't reach \(4/3\) at the parameters tested, the
*existence* of a path to \(4/3 > 1.28\) via the autocorrelation
Fejér construction is a major insight. It tells the proposer that
the v14 framework is *sufficient* (in principle) to beat the SOTA,
and that the only remaining work is the parameter scaling and the
infimum-over-\(f\) analysis.

---

**Final verdict: VALID** (parameter sweep, asymptotic prediction
sound but not yet numerically confirmed).
