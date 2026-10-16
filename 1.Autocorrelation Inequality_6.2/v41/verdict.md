# Verdict — Autocorrelation Inequality 6.2, v41

## 1. Summary

v41 finally executes the **\(K\)-sweep** that I have been
recommending in five consecutive verdicts (v33, v34, v35, v37, v39).
The two queued runs:

| \(N\) | \(K\) | \(n_{\text{tests}}\) | \(\Omega^{(41)}\) |
|---|---|---|---|
| 250 | 256 | 129 | (in progress) |
| 500 | 256 | 129 | (in progress) |

**No completed results yet.** The best rigorous bound on \(C_{6.2}\)
remains v37's \(1.2859\).

## 2. Verdict

**VALID-but-empty** (placeholder, but the experimental direction is
finally the right one).

## 3. Comments

### Comment 1. The proposer finally took the \(K\)-sweep recommendation.

This is the **right experiment**, queued at the right parameter
settings (\(K = 256\) — the most aggressive \(K\) value I
recommended). The proposer is also using a sensible
\(n_{\text{tests}} = 129\) (slightly more than \(K/2 = 128\), i.e.,
catching the Fejér mean's ~256 extrema with safety margin).

**Predicted outcome based on the v33–v37 \(K\)-sensitivity data:**

* \((250, 256)\): \(\sim 1.295\)–\(1.300\)
* \((500, 256)\): \(\sim 1.298\)–\(1.305\)

These are *substantially* better than the \(K = 64\) ceiling at
\(\sim 1.287\). If the predictions hold, the new bound will be
\(\sim 1.30\), nearly **doubling the SOTA exceedance margin**.

### Comment 2. Compute estimate.

At \(K = 256\), the SDP block size is \(\sim 514 \times 514\) (the
\(M\) matrix is \((2K+2)\times(2K+2)\)). MOSEK's per-iteration cost
is roughly \(O(K^{3})\), so doubling \(K\) from 64 to 256 (a 4×
increase) costs \(\sim 64\times\) more per iteration. But the LP
also has fewer \(N\)-dependent variables here, so the total compute
should be on the order of 30–60 minutes per run, not hours.

### Comment 3. Inheritance is clean.

v41 inherits from v30; no new bugs. The validity argument is
unchanged.

## 4. Feedback to the proposer

Please **wait for the runs to complete** and ship a single
complete v41 (or v42) with the numerical values.

Specifically:

1. **Report the completed numbers** with MOSEK status (optimal /
   feasible / time).
2. **Compare to the \(K = 64\) saturation curve** at \(N = 250, 500\).
3. **If the bound exceeds 1.29**, the project's \(K\)-sweep
   recommendation is vindicated, and the proposer should run
   \((1000, 256)\), \((250, 512)\), and \((100, 512)\) to push
   further.
4. **If the bound is similar to \(K = 64\) at the same \(N\)**,
   then the \(K\)-sweep also saturates near 1.287, and the
   v30 framework has its true ceiling there.

Either outcome is informative.

---

**Final verdict: VALID-but-empty** (placeholder for the long-awaited
\(K\)-sweep at \(K = 256\); best bound remains v37's
\(C_{6.2}\ge 1.2859\). Looking forward to the completed
numerical results).
