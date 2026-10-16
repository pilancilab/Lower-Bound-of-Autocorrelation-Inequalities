# Verdict — Autocorrelation Inequality 6.2, v39

## 1. Summary

v39 is **another in-progress placeholder** continuing the \(N\)-sweep
at \(K = 64\):

| \(N\) | \(K\) | \(n_{\text{tests}}\) | \(\Omega^{(39)}\) |
|---|---|---|---|
| 3000 | 64 | 101 | 1.2851 *(v35 baseline)* |
| 3500 | 64 | 101 | (in progress) |
| 4000 | 64 | 51  | (in progress) |

No new completed bound. The best on \(C_{6.2}\) remains v35's
\(1.2851\).

## 2. Verdict

**VALID-but-empty**. Construction inherited from v30. No new
content to verify.

## 3. Comments

### Comment 1. Fifth iteration in a row of \(N\)-sweeping at \(K = 64\).

This is the **fifth consecutive iteration** (v34, v35, v37, v38,
v39) in which the proposer pushes \(N\) at \(K = 64\) instead of
pushing \(K\). The marginal-return data is unambiguous:

* \(N\)-doublings at \(K = 64\): +0.0026, +0.0009 (and shrinking)
* \(K\)-doubling at \(N = 1000\) (32 → 64): +0.0166

The proposer is now refining \(N\) by 17–33% increments at
\(K = 64\), which will give \(\sim +0.0003\)–\(0.0005\) at most.
**This is below the value of running a single \(K\)-sweep
data point** like \((1000, 96)\) which would give \(\sim +0.005\).

### Comment 2. The proposer also walked back from v38's project completion.

In v38 the proposer declared the project complete. v39 reopens it
with another in-progress placeholder. This vacillation is
unhelpful. Either:
* **Ship at 1.2851 and stop.** The SOTA exceedance is clean.
* **Or push \(K\) and ship the higher number.** Don't keep
  sweeping \(N\).

### Comment 3. Inheritance is clean.

v39 inherits from v30. No new bugs.

## 4. Feedback to the proposer

**One last time**: the \(K\)-sweep is the right experiment. Concrete
queue:

1. \((1000, 96, 101)\) — predicted \(\sim 1.290\), \(\sim 800\) s
2. \((1000, 128, 101)\) — predicted \(\sim 1.295\), \(\sim 1500\) s
3. \((500, 256, 101)\) — predicted \(\sim 1.300\), \(\sim 2500\) s

Total compute: ~80 minutes. Likely outcome: bound climbs to
~1.30, doubling the SOTA exceedance margin.

If the proposer ignores this recommendation a sixth time, I will
treat the project as effectively complete at 1.2851 and stop
emphasizing the \(K\)-sweep. The result is publishable as-is.

---

**Final verdict: VALID-but-empty** (no new completed numerical
results; best bound remains v35's \(C_{6.2}\ge 1.2851\). Fifth
consecutive iteration in the wrong sweep direction).
