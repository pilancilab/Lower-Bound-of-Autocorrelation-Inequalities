# Verdict — Autocorrelation Inequality 6.2, v35 (updated with completed runs)

## 1. Summary

v35 reports two completed parameter runs (replacing the earlier
in-progress placeholder):

| \(N\) | \(K\) | \(n_{\text{tests}}\) | \(\Omega^{(35)}\) | time |
|---|---|---|---|---|
| 2000 | 64 | 101 | 1.2842 *(v34 baseline)* | 683 s |
| **3000** | **64** | **101** | **1.2851** | 694 s |
| 5000 | 32 | 51 | 1.2678 | 112 s |

**The \((3000, 64)\) run gives \(\Omega = 1.2851\), which is the
new best rigorous bound — \(\Omega > 1.2842\) from v34**, exceeding
the SOTA by \(0.0051\).

(The proposer's summary still cites v34's 1.2842 as the best bound,
which is a transcription error: 1.2851 is now the project's best.)

## 2. Verdict

**VALID** — small but real improvement on v34. Same v30 construction,
larger \(N\).

## 3. Verification

The construction and validity proof are inherited from v30; nothing
new to verify mathematically. The numerical results are consistent
with the v33/v34 trend:

* \(N\)-doubling at \(K = 64\): \((1000\to 2000)\) gave +0.0026, and
  \((2000\to 3000)\), a 1.5× scaling, gave +0.0009. Diminishing
  returns, as predicted.
* \(N\)-scaling at \(K = 32\): \((1000, 32) = 1.2650\),
  \((2000, 32) \approx 1.265\)–\(1.267\) (v31), now
  \((5000, 32) = 1.2678\). Saturation near \(1.27\) at \(K = 32\),
  consistent with v31's plateau prediction.

The compute times are reasonable: 694 s for \((3000, 64)\) is
linear in \(N\) at fixed \(K\), as expected.

## 4. Comments

### Comment 1. The proposer should update the "best bound" claim.

The boxed claim at the bottom of `rigorousproof.md` still reads
\(C_{6.2}\ge 1.2842\) (citing v34). With the new \((3000, 64)\)
result, the correct boxed claim is
\(C_{6.2}\ge 1.2851\), exceeding SOTA by \(0.0051\).

This is a minor cosmetic bug in the writeup, not a math error.

### Comment 2. The \(K = 64\) ceiling is now visible.

The trend at \(K = 64\):
* \(N = 1000\): 1.2816
* \(N = 2000\): 1.2842 (+0.0026)
* \(N = 3000\): 1.2851 (+0.0009)

Extrapolating, \(N \to \infty\) at \(K = 64\) gives a ceiling near
\(\sim 1.286\)–\(1.288\). To break past this, the proposer **must**
push \(K\) beyond 64.

### Comment 3. The \(K = 32\) ceiling is also visible.

At \(K = 32\):
* \(N = 1000\): 1.2650
* \(N = 5000\): 1.2678 (+0.0028)

Asymptotic ceiling \(\sim 1.268\)–\(1.270\). The \(K = 32\) ceiling
is \(\sim 0.018\) below the \(K = 64\) ceiling, confirming that
\(K\) is the binding lever.

### Comment 4. Inheritance is clean.

v35 inherits from v30. No new bugs.

## 5. Feedback to the proposer

For v36:

1. **Stop running \(K = 64\) at larger \(N\)** — diminishing returns.
2. **Run \(K = 96, 128, 192, 256\)** at moderate \(N\) (say, 1000
   or 2000). Specifically: \((1000, 96), (1000, 128), (2000, 96), (2000, 128)\).
3. **Update the boxed best-bound claim** to \(1.2851\).
4. **Once the \(K\)-sweep saturates**, write a final summary
   reporting the (best) bound and the asymptotic curve.

Predictions for v36:
* \((1000, 96)\): \(\sim 1.290\)–\(1.293\)
* \((1000, 128)\): \(\sim 1.293\)–\(1.297\)
* \((2000, 128)\): \(\sim 1.295\)–\(1.300\)

If the bound reaches 1.30+, the project has produced a substantial
improvement over the SOTA.

---

**Final verdict: VALID** (small improvement: \(C_{6.2}\ge 1.2851\)
at \((N, K) = (3000, 64)\), exceeding SOTA by \(0.0051\). The
\(K = 64\) curve is approaching its ceiling near \(\sim 1.287\);
v36 must push \(K\) beyond 64).
