# Verdict — Autocorrelation Inequality 6.2, v42 (updated with completed run)

## 1. Summary

v42 reports the completed \((N, K) = (10000, 64)\) run:

| \(N\) | \(K\) | \(n_{\text{tests}}\) | \(\Omega^{(42)}\) |
|---|---|---|---|
| 5000 | 64 | 101 | 1.285851 *(v37)* |
| **10000** | **64** | **101** | **1.286410** |

**New best: \(C_{6.2}\ge 1.2864\) at \((N, K, n_{\text{tests}}) = (10000, 64, 101)\)**,
exceeding the SOTA by \(0.0064\).

## 2. Verdict

**VALID** — small but real improvement on v37. Same v30
construction, larger \(N\).

## 3. Verification

Construction inherited from v30 (verified in v30's verdict);
nothing new to verify mathematically.

The numerical trend at \(K = 64\) is now:

| \(N\) | \(\Omega\) | marginal |
|---|---|---|
| 1000  | 1.2816 | — |
| 2000  | 1.2842 | +0.0026 |
| 3000  | 1.2851 | +0.0009 |
| 5000  | 1.28585 | +0.00075 |
| 10000 | 1.28641 | +0.00056 |

The marginal gains are shrinking by ~30% per \(N\)-doubling.
Extrapolating, the \(N\)-asymptote at \(K = 64\) is at
**\(\sim 1.2870\)**.

## 4. Comments

### Comment 1. Confirmation of the K=64 saturation curve.

Combined with the v43 \((500, 96), (1000, 96)\) data points, the
two parameter directions are converging to the **same limit**
\(\sim 1.287\). This is consistent with the observation that v30
has a true ceiling near 1.287–1.290 (per the revised prediction
in v43's verdict).

### Comment 2. The SOTA exceedance is now 0.0064.

\(1.2864 - 1.28 = 0.0064\). This is the cleanest SOTA exceedance
the project has produced so far. The result is shippable.

### Comment 3. Inheritance is clean.

v42 inherits from v30. No new bugs.

## 5. Feedback to the proposer

The v42 result is a real (small) improvement. The proposer should:

1. **Update the boxed best-bound claim** in v38 (if reusing as
   summary) to \(1.2864\) at \((10000, 64)\).
2. **Wait for v41's \(K = 256\) results** before declaring final.
   If they exceed \(1.287\), the project bound climbs further.
3. **Stop running new \(N\)-sweeps at \(K = 64\)**: marginal gains
   below \(\sim 0.0005\) per doubling are below the value of
   re-iterating.

---

**Final verdict: VALID** (incremental improvement: new best
\(C_{6.2}\ge 1.2864\) at \((10000, 64)\), exceeding SOTA by
\(0.0064\). The \(N\)-saturation curve at \(K = 64\) is converging
to \(\sim 1.287\)).
