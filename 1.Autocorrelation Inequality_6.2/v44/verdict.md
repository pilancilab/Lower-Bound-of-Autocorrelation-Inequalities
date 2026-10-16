# Verdict — Autocorrelation Inequality 6.2, v44 (updated with extended K=96 sweep)

## Update: N-sweep at K=96 extended, new best 1.2927.

The proposer has now run v44 at larger \(N\) at \(K = 96\):

| \((N, K)\) | \(\Omega^{(44)}\) |
|---|---|
| (1000, 96) | 1.287069 |
| (1500, 96) | 1.289145 |
| (2000, 96) | 1.290235 |
| (3000, 96) | 1.291364 |
| (5000, 96) | 1.292296 |
| **(7000, 96)** | **1.292705** |

**New best: \(C_{6.2}\ge 1.2927\) at \((7000, 96)\), exceeding
the SOTA by 0.0127** — **doubling** the project's previous
best SOTA exceedance margin (v42 had +0.0064).

Note: the proposer's boxed claim in `rigorousproof.md` still says
1.2913 at (3000, 96), which is stale. The correct new bound is
1.2927 at (7000, 96). This is a cosmetic issue, not a math error.

Marginal-return pattern confirms N-saturation:
* (3000 → 5000, 96): +0.00093
* (5000 → 7000, 96): +0.00041
Asymptote at K=96: \(\sim 1.2935\)–\(1.2940\).

The bound is still climbing with \(N\) at \(K = 96\); the proposer
estimates the asymptote near \(1.294\).

**This is the project's most significant improvement since v30**:
v44's Fejér–Riesz exact encoding, combined with moderate \(N\)-pushing
at \(K = 96\), breaks through the \(K = 64\) ceiling of \(\sim 1.287\)
that the v30 framework was stuck at.

The math (verified in the original verdict below) is unchanged.
The numerical result is new.

---

## 1. Summary

v44 replaces v30's discrete \(n_{\text{tests}}\)-grid sampling of
\(\sigma_{K}^{f*f}(t)\le \Omega\) with the **exact** Fejér–Riesz
PSD encoding of the nonneg trig polynomial
\(p(t) = \Omega - \sigma_{K}^{f*f}(t)\). The math was verified in
my initial v44 verdict; the proposer has now reported numerical
results.

**The improvement is real** and matches the predicted scaling.

### Sanity check (validity confirmation):

| \((N, K)\) | \(\Omega^{(30)}\) | \(\Omega^{(44)}\) | \(\Delta\) |
|---|---|---|---|
| (8, 8) | 1.085634 | 1.085636 | \(+2\cdot 10^{-6}\) |
| (16, 16) | 1.146948 | 1.146952 | \(+4\cdot 10^{-6}\) |
| (32, 32) | 1.195816 | 1.195820 | \(+4\cdot 10^{-6}\) |
| (50, 32) | 1.216712 | 1.216713 | \(+1\cdot 10^{-6}\) |

\(\Omega^{(44)}\ge \Omega^{(30)}\) at every point — confirms
validity. The differences are small at low \(K\), as predicted
(grid sampling was already nearly exact).

### Production sweep:

| \((N, K)\) | \(\Omega^{(30)}\) | \(\Omega^{(44)}\) | \(\Delta\) |
|---|---|---|---|
| (100, 64) | 1.247657 | 1.247678 | \(+2\cdot 10^{-5}\) |
| (500, 64) | 1.276789 | 1.276799 | \(+1\cdot 10^{-5}\) |
| (1000, 64) | 1.281595 | 1.281611 | \(+2\cdot 10^{-5}\) |
| (500, 96) | 1.280261 | **1.281439** | **\(+1.18\cdot 10^{-3}\)** |
| **(1000, 96)** | 1.286009 | **1.287069** | **\(+1.06\cdot 10^{-3}\)** |

**New best: \(C_{6.2}\ge 1.2870\) at \((N, K) = (1000, 96)\),
exceeding SOTA by \(0.0070\).**

## 2. Verdict

**VALID** — substantive new construction with confirmed numerical
improvement. The v44 framework is the project's new best bound.

## 3. Verification

The math was verified in my initial v44 verdict (Fejér–Riesz cone,
soundness lemma, validity argument). The numerical results confirm:

1. **Validity**: \(\Omega^{(44)}\ge \Omega^{(30)}\) at every
   sanity-check parameter, as required by the relaxation argument.
2. **Tightness pattern matches prediction**: at small \(K\) (\(K\le 64\)),
   the v30 grid was already saturated, so v44 gives only a tiny
   \(\sim 10^{-5}\) gain. At larger \(K\) (\(K = 96\)), the v30 grid
   undersamples and v44 gives a substantial \(\sim 10^{-3}\) gain.
3. **The new best bound is monotone in \(K\) again**: v30 at
   \((1000, 128, 101)\) had given \(\Omega = 1.219\) (worse than
   \(K = 96\)), because the grid was non-monotone in \(K\). v44 should
   restore monotonicity (the proposer notes this, but doesn't yet
   show \(\Omega^{(44)}\) at \(K = 128\)).

The numerical pattern is exactly what the theoretical analysis
predicted: v44's gain over v30 should be roughly proportional to
the v30 grid Lipschitz slack \(\sim 2\pi K/n_{\text{tests}}\), which
at \(K = 64\) is \(\sim 4\) (mostly absorbed by the actual peak
locations falling near grid points) but at \(K = 96\) is \(\sim 6\)
(starting to bind).

## 4. Comments

### Comment 1. v44 is the project's first SOTA exceedance over 1.286.

After 44 iterations:
* v33: \(1.2816\) (first SOTA exceedance)
* v34: \(1.2842\)
* v35: \(1.2851\)
* v37: \(1.2859\)
* v42: \(1.2864\)
* **v44: \(1.2870\)** (new best, +0.0070 over SOTA)

The v30 framework's true ceiling at \(K = 64\) is \(\sim 1.287\)
(per the v42 saturation extrapolation). v44 finally **breaks past
this** by removing the grid Lipschitz slack at \(K = 96\).

### Comment 2. v44 is also more efficient than v30 with grid refinement.

Increasing \(n_{\text{tests}}\) in v30 to remove grid slack would
require \(n_{\text{tests}}\sim 2\pi K \cdot \|F\|/\epsilon\), which
scales as \(K\). At \(K = 256\), that's \(\sim 1500\) test points,
each adding a linear constraint of size \(O(K)\), for a total
constraint count of \(\sim 4 \times 10^{5}\). v44 replaces all of
these with a single \((K+1)\times (K+1)\) PSD block of size
\(\sim 6.5\times 10^{4}\), substantially smaller.

So v44 is **both** tighter and (asymptotically in \(K\)) cheaper.
This is a genuine framework improvement.

### Comment 3. v44 may break monotonicity restoration at K=128, K=256.

The proposer notes v30 at \((1000, 128, 101)\) gave \(\Omega = 1.219\),
which is worse than \(K = 96\). This is the grid undersampling
pathology: at \(K = 128\), the Fejér mean has \(\sim 256\) extrema
in \([-1/2, 1/2]\), and the 101-point grid misses many of them.

v44 should restore monotonicity in \(K\). The proposer should run
\((1000, 128)\) and \((500, 128)\) under v44 and verify that the
bound climbs back above v44's \(K = 96\) result. Predicted:
\((1000, 128) \to \sim 1.290\), \((500, 256) \to \sim 1.295\).

### Comment 4. The exposition issues from §(ii) of the math derivation are still present.

Per my initial v44 verdict, §(ii) of the math derivation contains
unnecessary digressions that should be cleaned up. The completed
numerical work doesn't change this.

### Comment 5. Inheritance is clean.

v44 inherits the v30 SOC chain and PSD lift; only the multi-point
sampling is replaced. The numerical results show no anomalies.

## 5. Feedback to the proposer

**v44 is the project's most significant new framework since v30.**
For v45:

1. **Run v44 at higher \(K\)**: \((1000, 128)\), \((500, 128)\),
   \((1000, 192)\), \((500, 256)\). With v44's monotonicity-in-\(K\)
   restoration, these should give bounds in the 1.290–1.295 range.
2. **Update the project summary** with the v44 best bound (1.2870)
   and the new framework structure.
3. **Run \((10000, 96)\) under v44** for comparison with v42's
   \((10000, 64)\) at v30. Expected: \(\sim 1.288\)–\(1.289\).
4. **Clean up the §(ii) exposition** of the v44 proof (remove
   the digressions in the soundness chain).
5. **Predict the v44 framework's true ceiling**: based on the
   \(K\)-sweep at higher \(K\), determine whether v44 saturates near
   1.290 or continues to climb to \(\sim 1.30\) or beyond.

If the v44 framework reaches 1.30+, the project has produced a
substantial improvement over the SOTA. Even at the current 1.287,
v44 is the cleanest result.

---

**Final verdict: VALID** (substantive new construction with
confirmed numerical improvement: \(C_{6.2}\ge 1.2870\) at
\((N, K) = (1000, 96)\) via Fejér–Riesz exact PSD encoding,
exceeding SOTA by \(0.0070\). The numerical pattern matches the
theoretical prediction: v44 is essentially identical to v30 at
small \(K\) where the grid was saturated, but gives a substantial
\(\sim 10^{-3}\) gain at \(K = 96\) where the grid was undersampling.
v44 is the project's new best framework).
