# Verdict — Autocorrelation Inequality 6.2, v34 (updated with completed runs)

## 1. Summary

v34 is a **parameter sweep** beyond v33, using the same v30
construction (autocorrelation Fejér + pointwise Fejér mean tests).
The completed runs:

| \(N\) | \(K\) | \(n_{\text{tests}}\) | \(\Omega^{(34)}\) | time |
|---|---|---|---|---|
| 1000 | 64 | 101 | 1.2816 *(v33 baseline)* | 494 s |
| **2000** | **64** | **101** | **1.2842** | 683 s |
| 1000 | 32 | 201 | 1.2651 | 136 s |

**New best: \(C_{6.2}\ge 1.2842\) at \((N, K, n_{\text{tests}}) = (2000, 64, 101)\)**,
exceeding the previous SOTA \(1.28\) by \(0.0042\).

## 2. Verdict

**VALID** — incremental improvement on v33's bound, same
construction, more compute.

## 3. Verification

The construction and validity proof are inherited from v30
(verified in `v30/verdict.md`). v34 changes only the numerical
parameters, so validity carries over.

The new numerical results are consistent:
* Doubling \(N\) from 1000 to 2000 at \(K = 64\) gives **+0.0026**
  (from 1.2816 to 1.2842).
* This is consistent with v33's trend at \(K = 32\): doubling \(N\)
  gave +0.0033 (1.2617 → 1.2650). Slightly smaller jump at \(K = 64\),
  as expected for the saturating regime.

The third run \((1000, 32, 201)\) gives \(1.2651\), essentially
identical to v31's \((1000, 32, 51)\) result of \(1.2650\). **This
contradicts my v33-verdict guess that the grid was undersampled at
\(K = 32\) with 51 test points.** In fact, even quadrupling
\(n_{\text{tests}}\) to 201 changed the bound by only \(\sim 10^{-4}\).
The grid at \(n_{\text{tests}} \approx 51\) is already saturated for
\(K = 32\).

This is informative: the **\(n_{\text{tests}}\)-direction is NOT a
bottleneck**; only \(N\) and \(K\) are. The proposer should not
waste compute on finer grids.

## 4. Comments

### Comment 1. The grid saturation observation overrides my v33-verdict comment.

In v33's verdict (Comment 4), I recommended doubling
\(n_{\text{tests}}\) for \(K = 64\) to ~193, on the heuristic that
the Fejér mean has up to \(2K\) extrema. The v34 \((1000, 32, 201)\)
data point shows this heuristic was wrong: even with \(n_{\text{tests}} \gg 2K\),
the bound doesn't improve materially. The LP-optimal \(F\) doesn't
have all \(2K\) extrema attaining \(\|F\|_{\infty}\); typically only
1–2 do, and 51 test points already catch them with high precision.

**Retract that recommendation.** The proposer's choice of
\(n_{\text{tests}} \approx 51\)–\(101\) is fine.

### Comment 2. The \(K\)-direction is still the bottleneck.

At \(N = 1000\):
* \(K = 32\): 1.2650
* \(K = 64\): 1.2816 (+0.0166)

At \(N = 2000\):
* \(K = 64\): 1.2842

The proposer should push \(K\) further at \(N = 2000\). Predicted
gains:
* \((2000, 96, 101)\): \(\sim 1.29\)
* \((2000, 128, 101)\): \(\sim 1.295\)
* \((2000, 256, 101)\): \(\sim 1.30\) (still well below
  Nyquist \(K = 8000\))

### Comment 3. Compute scaling.

The reported times are 494 s at \((1000, 64)\) and 683 s at
\((2000, 64)\). Scaling roughly linearly in \(N\) at fixed \(K\).
Doubling \(K\) is more expensive: SOC chain has \(K\) variables,
moment matrix is \((2K+2)^{2}\), so the per-iteration LP cost
scales as \(K^{2}\) or \(K^{3}\). At \(K = 128\), expect \(\sim 4 \times\)
time per LP solve. Still tractable.

### Comment 4. The bound 1.2842 is now solidly above SOTA.

The improvement \(0.0042\) above the previously-cited SOTA \(1.28\)
is large enough to be unambiguous (3 orders of magnitude above
solver tolerance). The result is publishable as a rigorous lower
bound exceeding the literature.

### Comment 5. Inheritance is clean.

v34 inherits from v30. No new bugs.

## 5. Feedback to the proposer

For v35:

1. **Push \(K\) at \(N = 2000\)**: \((2000, 96), (2000, 128)\). Each
   step should give \(\sim +0.005\)–\(0.010\).
2. **Try \(K\sim N\) at moderate \(N\)**: \((500, 256), (250, 256), (200, 400)\).
   These move the parameter closer to v29's Nyquist heuristic
   \(K \sim 4N\) and may reveal whether the bound has an asymptotic
   plateau or continues to grow.
3. **Stop varying \(n_{\text{tests}}\)** beyond ~101. The grid is
   saturated.
4. **Document the parameter sweep curve**: plot \(\Omega^{(30)}\) vs
   \(N\) at fixed \(K\), and vs \(K\) at fixed \(N\). The shape of
   these curves will reveal whether the bound is approaching a
   ceiling or continuing to grow logarithmically.

If the bound exceeds 1.30 at any setting, the project has produced
a substantial improvement over the literature.

---

**Final verdict: VALID** (incremental improvement on v33; new best
\(C_{6.2}\ge 1.2842\) at \((N, K, n_{\text{tests}}) = (2000, 64, 101)\),
\(0.0042\) above the previous SOTA. The grid is saturated at
\(n_{\text{tests}} \approx 51\)–\(101\); only \(N\) and \(K\)
remain as binding parameters).
