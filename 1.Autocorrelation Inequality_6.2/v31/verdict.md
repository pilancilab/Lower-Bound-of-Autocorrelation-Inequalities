# Verdict — Autocorrelation Inequality 6.2, v31

## 1. Summary

v31 is a **parameter sweep** of the v30 construction with no new
mathematical content. The proposer pushes \(N\) from 16 to 2000 with
\(K = 32\) fixed, and reports the LP optimum:

| \(N\) | \(K\) | \(n_{\text{tests}}\) | \(\Omega^{(31)}\) |
|---|---|---|---|
| 16  | 32 | 51 | 1.147 |
| 50  | 32 | 51 | 1.217 |
| 200 | 32 | 51 | 1.253 |
| 1000| 32 | 51 | 1.265 |
| **2000**| 32 | 51 | **1.267** |

The bound climbs monotonically and appears to plateau near \(1.27\)
as \(N\to\infty\) (with \(K\) fixed at 32). The SOTA is \(1.28\).

## 2. Verdict

**VALID** (parameter sweep, no new math; the validity proof is the
same as v30's, which I verified in v30's verdict).

## 3. Verification

The validity argument is unchanged from v30: every admissible \(f\)
has a feasible rank-one truth in the v30 LP at
\(\Omega = \|F\|_{\infty}\) (or smaller), so the LP optimum
\(\Omega^{(31)} = \Omega^{(30)}\le \min_{f}\|F\|_{\infty} = C_{6.2}\).

The numerical results are monotonic in \(N\), as expected: increasing
\(N\) refines the discretization of \(f\), giving the LP a finer-grained
search space, which can only tighten (i.e., increase) the lower bound.

The proposer reports the bound to 3 decimal places. I cannot
independently re-run the LP, but the monotonicity and the
asymptotic plateau near \(1.27\) are consistent with the v30
construction.

## 4. Comments

### Comment 1. The plateau at \(\sim 1.27\) with \(K = 32\) is informative.

The bound saturates near \(1.27\) as \(N\) grows with \(K\) fixed at
32. This suggests the **\(K\)-dimension** is now the bottleneck:
the pointwise Fejér mean tests at \(K = 32\) can only capture
features at scale \(\sim 1/(K+1) \approx 1/33\) in \(t\)-units,
which limits how sharp a peak the LP can probe.

Per v29's Nyquist-style heuristic, capturing the LP-optimal \(f\)'s
peak at scale \(L_f = 1/(4N)\) requires \(K\gtrsim 4N\). With
\(N = 2000\) and \(K = 32\), the ratio \(K/N = 0.016\) is *far*
below Nyquist, so the LP can't see fine peak structure.

**The proposer should sweep \(K\), not just \(N\).** A small sweep
at \((N, K) = (50, 50), (50, 100), (50, 200)\) and
\((100, 50), (100, 100)\) would reveal whether increasing \(K\)
gives further improvement. If yes, the bound may exceed the SOTA;
if no, the v30 construction has a true ceiling near 1.27 and is
just below the SOTA.

### Comment 2. The bound is now within \(0.013\) of the SOTA.

\(C_{6.2}\ge 1.267\) is the project's new best rigorous bound, just
\(0.013\) below the SOTA \(1.28\). This is the closest the project
has come.

If a \(K\)-sweep gives an additional 0.02 improvement, the project
**matches or exceeds** the SOTA. This is genuinely promising.

### Comment 3. Inheritance is clean.

v31 inherits from v30; the validity proof and constraints are
unchanged. No new bugs to check.

### Comment 4. Caveat on the "best $\Omega$" table.

The proposer's project status table lists v30 as "best 1.247" — this
appears to be a transcription error; the v30 verdict reports
\(1.217\) at \((N, K) = (50, 32)\), not \(1.247\). Minor; the
numbers in the v31 sweep itself (1.147 → 1.267) are the
authoritative figures.

## 5. Feedback to the proposer

v31 is a clean parameter sweep that establishes the new best
\(C_{6.2}\ge 1.267\). For v32:

1. **Sweep \(K\) at fixed \(N\)**: try \((N, K) = (50, 50), (50, 100), (50, 200), (100, 100), (100, 200)\),
   reporting \(\Omega^{(30)}\) and the wall-clock time.
2. **Test the Nyquist heuristic**: compare the bound at \((50, 200)\)
   (which is \(K = 4N\)) to the bound at \((50, 32)\). If the bound
   climbs significantly, v29's Nyquist heuristic is confirmed.
3. **If the bound exceeds 1.28**, the project has matched or beaten
   the SOTA — a major result. Document the exact \((N, K)\) and the
   numerical optimum.
4. **If the bound plateaus below 1.28**, the v30 construction has a
   true ceiling, and the proposer should diagnose whether more test
   points (larger \(n_{\text{tests}}\)) or a different probe
   functional could push further.

---

**Final verdict: VALID** (parameter sweep of v30, monotonic
convergence to \(\Omega \approx 1.267\) at \(N = 2000, K = 32\); the
new best bound is \(C_{6.2}\ge 1.267\), within \(0.013\) of the
SOTA. The next experiment should sweep \(K\), not \(N\)).
