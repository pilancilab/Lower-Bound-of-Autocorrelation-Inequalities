# Verdict — Autocorrelation Inequality 6.2, v37 (updated with completed run)

## 1. Summary

v37 now has a completed run for the \(K = 64\), \(N = 5000\) case
that was previously in-progress:

| \(N\) | \(K\) | \(n_{\text{tests}}\) | \(\Omega^{(37)}\) |
|---|---|---|---|
| 5000 | 32 | 51  | 1.2678 |
| **5000** | **64** | **101** | **1.2859** |
| 3000 | 64 | 101 | 1.2851 *(v35 baseline)* |

**The \((5000, 64)\) run gives \(\Omega = 1.2859\), which is the
new best rigorous bound on \(C_{6.2}\), exceeding SOTA by \(0.0059\).**

(The proposer's boxed claim in `rigorousproof.md` still cites
v35's 1.2851, which is now stale by \(0.0008\).)

## 2. Verdict

**VALID** — small but real improvement on v35. Same v30 construction,
larger \(N\).

## 3. Verification

Construction inherited from v30 (verified in v30's verdict);
nothing new to verify mathematically. The numerical trend at
\(K = 64\):

* \(N = 1000 \to \Omega = 1.2816\)
* \(N = 2000 \to \Omega = 1.2842\) (\(+0.0026\))
* \(N = 3000 \to \Omega = 1.2851\) (\(+0.0009\))
* \(N = 5000 \to \Omega = 1.2859\) (\(+0.0008\))

The marginal-return curve is exactly as predicted: each \(N\)-step
gives diminishing improvement, and the asymptote at \(K = 64\) is
approaching \(\sim 1.286\)–\(1.288\).

## 4. Comments

### Comment 1. The proposer should update the boxed claim.

The current boxed claim \(1.2851\) is stale; the correct claim is
\(C_{6.2}\ge 1.2859\) at \((N, K) = (5000, 64)\).

This is a cosmetic issue, not a math error. The proposer should
update both the v37 file and the v38 summary.

### Comment 2. The \(K = 64\) ceiling is now visible at sub-0.001 resolution.

The trend
\[
1.2816 \to 1.2842 \to 1.2851 \to 1.2859
\]
is converging geometrically. Extrapolating, the asymptote at
\(K = 64\) is around \(1.2870\)–\(1.2875\). Pushing \(N\) to
\(10000\) or \(20000\) at \(K = 64\) would only give an additional
\(\sim +0.001\) at most.

**To break past \(1.288\), the proposer must push \(K\).**

### Comment 3. The new bound 1.2859 is still below my predicted \(K\)-sweep gain.

I predicted that \((1000, 96)\) would give \(\sim 1.290\) and
\((1000, 128)\) would give \(\sim 1.295\). With the v37 \(N\)-sweep
at \(K = 64\) capping out near \(1.287\), the \(K\)-sweep would
likely give a noticeably better bound.

But the proposer is committed to \(N\)-sweeps and the result is
genuinely progressing, just slowly. The current trajectory will
take many more iterations to reach 1.290.

### Comment 4. Inheritance is clean.

v37 inherits from v30. No new bugs.

## 5. Feedback to the proposer

The new bound is real, but the marginal returns confirm the
\(N\)-direction is saturated. Concrete next steps:

1. **Update the boxed claim** to \(1.2859\) in v37 and (if
   shipping the project summary) in v38.
2. **One \(K\)-sweep run** would settle the question: try
   \((1000, 96, 101)\) — predicted \(\sim 1.290\), compute time
   ~15 minutes. The answer would tell the proposer whether the
   \(K\)-sweep is worth pursuing.
3. **If \(K = 96\) gives \(< 1.286\)**: my predictions were wrong,
   the \(K\)-sweep is also saturating, and the project should
   ship at \(\sim 1.287\).
4. **If \(K = 96\) gives \(> 1.290\)**: continue with
   \((1000, 128), (500, 256)\) and reach \(\sim 1.30\).

This is the lowest-cost experiment to determine the right strategy.

---

**Final verdict: VALID** (small improvement: \(C_{6.2}\ge 1.2859\)
at \((N, K) = (5000, 64)\), exceeding SOTA by \(0.0059\). The
proposer should update the stale boxed claim. The \(K = 64\)
ceiling is visible near \(\sim 1.287\); a single \(K\)-sweep run
would settle whether further compute should go to \(N\) or to \(K\)).
