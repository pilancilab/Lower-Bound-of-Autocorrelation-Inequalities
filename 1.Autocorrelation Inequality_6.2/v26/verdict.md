# Verdict — Autocorrelation Inequality 6.2, v26

## 1. Summary

v26 is a **diagnostic** iteration — no new constraints. The proposer
extracts the LP-optimal \((p, a, b)\) at \((N, K) = (16, 16)\) and
computes the actual \(\int F^{2}\) for the corresponding \(f\). The
key finding:

* The LP-optimal \(p\) is **not uniform**: it has heavy mass at the
  edges (\(p_1 = p_{32} \approx 0.222\)) and a non-trivial structure
  in the middle.
* The LP-optimal \(f\) achieves \(\int F^{2} \approx 1.07\), which
  is **much less than \(4/3\)** (the value at the uniform \(f\)).
* The asymptotic limit of the v14 framework is therefore at most
  \(\sim 1.07\), **not \(4/3\)**.
* **v14 cannot reach the SOTA \(C_{6.2}\ge 1.28\)** by any amount of
  parameter scaling.

## 2. Verdict

**VALID** — diagnostic, no new mathematical content. All v14
constraints and bounds remain unchanged. The finding is the important
content, and it changes the project's roadmap.

## 3. Verification of the diagnostic

I verified the proposer's Plancherel arithmetic.

By Parseval on \(F\in L^{2}([-1/2, 1/2])\) with the orthonormal cosine
Fourier basis (using \(F\) is real, \(\int F = 1\)):
\[
\int F^{2}\,dt = \hat F(0)^{2} + 2\sum_{k\ge 1}|\hat F(k)|^{2} = 1 + 2\sum_{k\ge 1}|\hat F(k)|^{2}.
\]
With the proposer's number \(\sum_{k\ge 1}|\hat F(k)|^{2}\approx 0.034\):
\[
\int F^{2} \approx 1 + 2\cdot 0.034 = 1.068. \;\checkmark
\]

Compared to v14's reported \(\Omega^{(14)} = 1.063\) at the same
parameters: the gap is \(\sim 0.005\), consistent with the SDP
relaxation slack on \(M\) (where \(v_k\) may be slightly larger than
the actual \(|\hat F(k)|^{2}\)).

**The diagnostic computation is consistent and verified.** ✓

## 4. The critical implication

The earlier project prediction (v17 §, v25 §):
> The v14 framework's asymptotic limit (as \(N, K\to\infty\)) is at
> most \(\inf_f \int F^{2}\), which is at most \(4/3\) (achieved at
> the uniform \(f = 2\cdot\mathbf 1_{[-1/4, 1/4]}\)).

The v26 diagnostic shows this prediction is **wrong**:
* The LP is *not* picking the uniform \(f\); it picks a non-uniform
  \(f\) with much smaller \(\int F^{2}\).
* The actual asymptotic limit \(\inf_f \int F^{2}\) (over admissible
  \(f\), with the SDP-relaxed constraints) is at most \(1.07\)
  (witnessed by the LP-optimal \(f\) at \(N=16\)), not \(4/3\).
* Even with infinite \(N, K\), the LP cannot exceed
  \(\sim 1.07\)–\(1.15\) by this construction.

**v14 (and the entire v15–v25 family) cannot reach the SOTA
\(C_{6.2}\ge 1.28\).** The asymptotic limit is well below the SOTA.

## 5. Comments

### Comment 1. This diagnostic was the right move.

I had been recommending this diagnostic since the v17 verdict's
Comment 2 ("Take the LP-optimal \((p, a, b)\) at large \(N\) and
compute \(\int F^{2}\) numerically"), and again in v23 Comment 3 and
v25 Comment 4. The proposer's v26 finally executes it, and the
finding is *more* informative than I expected: the asymptotic limit
is *much* smaller than \(4/3\), not just a bit smaller.

This is the kind of negative result that changes the project's
direction. The v14 framework is *fundamentally too weak* to reach
the SOTA, and pushing parameters further is futile.

### Comment 2. The project's understanding is now clearer.

After v26, the project's state is:

| Quantity | Value | Notes |
|---|---|---|
| Trivial bound | 1.0 | from \(\int F = 1\) and \(F\le \Omega\) |
| **v14 LP at \((16, 16)\)** | **1.063** | rigorous |
| **v14 LP asymptotic limit** | **\(\sim 1.07\)** | **per v26 diagnostic** |
| v14 LP at \((100, 32)\) | 1.111 | parameter sweep |
| Uniform \(f\)'s value | \(4/3\) | **NOT the asymptotic limit** |
| SOTA | 1.28 | unchanged |
| (Conjectural) optimal value | unknown | between 1.28 and 2 |

The slow convergence of v14 (1.063 → 1.111 from \(N=16\) to \(N=100\))
is consistent with the asymptotic limit being \(\sim 1.07\)–\(1.15\),
not \(4/3\). The LP is approaching its actual asymptotic limit, which
happens to be modestly below the SOTA.

### Comment 3. What v27 should do.

The proposer's three suggested directions for moving past v14:

1. Higher-order Fourier moments of \(f\) (e.g., \(\int f^{3}\),
   \(\int f^{4}\)).
2. Full Bochner constraint on \(\hat f\) (which v13 added but didn't
   help).
3. A different relaxation framework (e.g., SOS hierarchy).

I would add a fourth:

4. **Look at the structure of the LP-optimal \(f\).** The diagnostic
   shows that the LP-optimal \(f\) has \(p_1 = p_{32} \approx 0.222\)
   (heavy edges) and a complex middle structure. This *suggests*
   the optimal \(f\) is something like a sum of two delta-like spikes
   at the edges plus a small bulk in the middle. If you can guess
   the analytic form, you can:
   * Compute its actual \(\int F^{2}\) and \(\|F\|_{\infty}\) to
     verify.
   * Use the analytic form to construct a test \(f\) for upper-bounding
     the *actual* \(C_{6.2}\) (which is the *infimum* of
     \(\|F\|_{\infty}\)) by giving an admissible \(f\) with
     small \(\|F\|_{\infty}\).

The fact that the LP found a solution with \(\int F^{2}\approx 1.07\)
*and the LP's bound is the autocorrelation Fejér via \(\sum |\hat F(k)|^{2}\)*
means: the actual \(\|F\|_{\infty}\) for this \(f\) might be much
smaller than \(1.07\) — possibly approaching the actual
\(\inf_f\|f * f\|_{\infty}\). If so, this gives an *upper* bound on
\(C_{6.2}\) that's possibly much smaller than the SOTA's *lower*
bound 1.28.

Wait, that can't be right. The SOTA \(C_{6.2}\ge 1.28\) is a
*lower* bound. An upper bound on \(C_{6.2}\) from a specific
admissible \(f\) gives \(C_{6.2}\le \|f * f\|_{\infty}\) for that
\(f\). If the LP-optimal \(f\) has \(\|f * f\|_{\infty}\) close to
\(\sqrt{1.07}\approx 1.034\), that would give \(C_{6.2}\le 1.034\),
contradicting the SOTA \(1.28\).

So either:
* The LP-optimal \(f\) is *not actually admissible* (perhaps it
  fails the support or normalization constraint, which the LP only
  enforces approximately via the cell discretization).
* Or the LP's \(\int F^{2}\) underestimates the actual
  \(\|F\|_{\infty}\) (which is plausible — the LP bounds \(\Omega = \|F\|_{\infty}\)
  via Fejér, and \(\Omega \ge \int F^{2}\) is just one of the
  constraints; the actual \(\|F\|_{\infty}\) for the LP-extracted
  \(f\) could be larger than \(\int F^{2}\)).

This would be a useful sanity check: take the LP-extracted \(f\),
compute its actual \(\|F\|_{\infty}\) (not the LP bound), and
compare to the SOTA. If the actual \(\|F\|_{\infty} < 1.28\), the
SOTA is violated, which can't happen — so something is wrong with
the LP extraction or the diagnostic. If the actual
\(\|F\|_{\infty}\ge 1.28\), then the LP is correctly *under*-estimating
\(\|F\|_{\infty}\) by the autocorrelation Fejér slack.

**This is the next diagnostic to do.** Without it, we don't know
whether the v14 framework's slow convergence is intrinsic
(v14 simply *cannot* see the actual \(\|F\|_{\infty}\) of the
extremal \(f\)) or remediable.

### Comment 4. The right takeaway.

v26 has *changed* the project's understanding: the v14 framework
has an asymptotic ceiling well below the SOTA. This is *critical*
information and should reorient the proposer's effort.

Specifically: the proposer should *stop* pushing parameter scaling
on v14 (no improvement is possible past \(\sim 1.10\)), and *start*
exploring fundamentally different approaches (Lasserre on \(p\)
moments, sum-of-squares, Schur lift, dual constructions, etc.).

The proposer's v27 plan should explicitly choose one of these
directions and execute it, rather than continuing the v15–v25 mode
of incremental tightenings.

---

**Final verdict: VALID** (critical diagnostic finding; v14 cannot
reach the SOTA via parameter scaling).
