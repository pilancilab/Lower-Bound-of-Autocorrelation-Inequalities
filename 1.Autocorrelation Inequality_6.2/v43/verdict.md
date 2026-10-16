# Verdict — Autocorrelation Inequality 6.2, v43 (updated with completed runs)

## 1. Summary

v43 now reports two completed \(K = 96\) runs:

| \(N\) | \(K\) | \(n_{\text{tests}}\) | \(\Omega^{(43)}\) |
|---|---|---|---|
| 500  | 96 | 101 | 1.2803 |
| **1000** | **96** | **101** | **1.2860** |

The proposer correctly notes that v42's
\((10000, 64) \to 1.2864\) is the overall best, slightly above
v43's \((1000, 96) \to 1.2860\). The two tracks (large \(N\)
at \(K = 64\) vs moderate \(N\) at \(K = 96\)) converge to the
**same limit near \(\sim 1.287\)**.

## 2. Verdict

**VALID** — informative \(K\)-sweep data; both \(K = 96\) bounds
are rigorous (v30 construction inherited).

## 3. Verification

Inherited validity from v30. The new data points refine the
picture from my v43-initial verdict.

### Cross-track comparison.

| Track A: \((N, 64)\) | \(\Omega\) |     | Track B: \((N, 96)\) | \(\Omega\) |
|---|---|---|---|---|
| \((500, 64)\) | 1.2768 | | \((500, 96)\) | 1.2803 |
| \((1000, 64)\) | 1.2816 | | \((1000, 96)\) | **1.2860** |
| \((2000, 64)\) | 1.2842 | | | |
| \((5000, 64)\) | 1.28585 | | | |
| \((10000, 64)\) | 1.28641 | | | |

The \(K = 96\) track at \(N = 1000\) gives \(1.2860\), which is
within \(0.0004\) of the \(K = 64\) track at \(N = 10000\)
(\(1.2864\)). **The two tracks are converging to the same limit
\(\sim 1.287\)**, as the v43-initial verdict predicted.

### Cost comparison.

* \((1000, 96)\): smaller \(N\), larger \(K\). The LP block size
  for \(M\) is \((2K+2)^{2}\), so \(K\)-cost is \(O(K^{2})\). At
  \(K = 96\) vs \(K = 64\), that's \(\sim 2.25\times\) the per-LP
  size. With \(N\) reduced by 10×, the total compute is comparable
  or less.
* So **\((1000, 96)\) is a more efficient way to reach \(\sim 1.286\)**
  than \((10000, 64)\). This vindicates v33's original \(K\)-sweep
  recommendation, just at a smaller magnitude than I initially
  predicted.

## 4. Comments

### Comment 1. The v43 result confirms the framework's true ceiling.

The two-track convergence to \(\sim 1.287\) is strong evidence that
v30's framework has a true asymptotic ceiling near
\(1.287\)–\(1.290\). To break past this, the proposer needs either:

* **Much larger \(K\)** (e.g., \(K = 256\), still in v41), or
* **A different relaxation framework** (e.g., higher-order
  Lasserre on \(f\)).

The v41 \(K = 256\) runs are the critical experiment.

### Comment 2. The proposer's framing is correct.

The proposer correctly notes that v42's bound \(1.2864\) is the
overall best, not v43's \(1.2860\). Good attention to detail.

### Comment 3. Inheritance is clean.

v43 inherits from v30; no new bugs.

## 5. Feedback to the proposer

The v43 + v42 combination is the project's strongest two-track
verification of the v30 framework's saturation curve. For v44+:

1. **Wait for v41 \(K = 256\) runs.** These are the only remaining
   plausibly-tight \(K\)-sweep experiments.
2. **If v41 gives bounds in the 1.290+ range**, the framework's
   ceiling is higher than \(1.287\); push further.
3. **If v41 gives \(\le 1.288\)**, the framework's true ceiling is
   indeed \(\sim 1.288\) and the project should ship at v42's
   \(1.2864\).
4. **Update v38's boxed claim** to \(1.2864\).

---

**Final verdict: VALID** (informative \(K\)-sweep data point at
\((1000, 96) = 1.2860\), confirming the two-track convergence to
\(\sim 1.287\). Best bound on \(C_{6.2}\) remains v42's \(1.2864\)
at \((10000, 64)\)).
