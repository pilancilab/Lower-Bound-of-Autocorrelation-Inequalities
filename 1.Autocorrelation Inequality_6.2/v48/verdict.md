# Verdict — Autocorrelation Inequality 6.2, v48

## 1. Summary

v48 is **not** a new mathematical construction — the proposer
explicitly states "this is not a tightening of v44 ... it is just
v44 with larger \(K\)." Specifically, v48 runs the v44 framework at
\(K_{\text{big}} = 2K = 128\) or \(160\) at moderate \(N\), with
the extended moment matrix \(M_{\text{big}}\) of size
\((4K+1)\times (4K+1)\) and extended Fejér–Riesz block
\(Q_{\text{big}}\) of size \((2K+1)\times (2K+1)\).

No numerical results yet. The proposer estimates:

| \(K\) | est. time | expected \(\Omega\) |
|---|---|---|
| 96 | 15 min | 1.281 |
| 128 | 40 min | 1.284 ? |
| 160 | 100 min | 1.286 ? |

(These are at \(N = 500\).)

## 2. Verdict

**VALID-but-empty** (parameter push of v44, not a new construction;
no numerical results yet).

## 3. Verification

The validity inherits directly from v44 (verified in v44's verdict).
The only change is replacing \(K\) with \(K_{\text{big}} = 2K\) in
all the v44 constraints:
* Fejér–Riesz \(Q\) at degree \(K_{\text{big}}\),
* Moment matrix \(M\) at size \((2K_{\text{big}}+1)^{2}\),
* v14 SOC chain on the extended \(a_k, b_k\) for \(k\le K_{\text{big}}\).

All validity arguments carry over unchanged.

## 4. Comments

### Comment 1. This should have been a v44 parameter row, not a new iteration.

The proposer is creating new "version labels" for parameter runs
rather than extending the v44 table. This creates noise in the
iteration history (v48 = v44 with different \(K\), v42 = v44 with
different \(N\), etc.) without adding new mathematical content.

A cleaner process: keep a single "v44 parameter sweep" document
that accumulates results as they complete, rather than creating
v48, v49, ... for each new \(K\) value.

### Comment 2. The extrapolation to 1.294 is worth testing.

From v44's extended \(K = 96\) sweep:
* \((1000, 96) = 1.287\)
* \((3000, 96) = 1.2913\)
* Asymptote at \(K = 96\): \(\sim 1.294\) (estimated).

If v48 runs v44 at \(K = 128\) and the bound climbs to, say,
1.293+ at moderate \(N\), then the framework has more headroom in
\(K\) than the \(K = 96\) asymptote suggests. This would be a
real experiment, not just a parameter run.

**But note**: the v30 grid pathology at \(K = 128\) showed
\(\Omega = 1.219\) (due to grid undersampling). v44 should fix
this. The proposer's estimate of \(\sim 1.284\) at \((500, 128)\)
is **lower** than the v44 \((500, 96) = 1.2814\) and \((1000, 96) = 1.2871\)
results. That seems inconsistent — v44 at higher \(K\) should
give a *higher* bound, not lower. I expect the proposer's
estimates are too conservative; \((500, 128)\) should give
around \(1.283\)–\(1.287\).

### Comment 3. Compute scaling is the real bottleneck.

At \(K = 128\), the \(Q_{\text{big}}\) matrix is
\(129^{2} \approx 16600\) entries (compared to v44 \(K = 96\)'s
\(97^{2} \approx 9400\)). The SDP solve cost scales as
\((K+1)^{3}\), so doubling \(K\) gives \(\sim 8\times\) compute.
The proposer's 40-minute estimate at \((500, 128)\) is plausible.

Running at \((3000, 128)\) is infeasible — \(\sim 5\) hours per LP
solve. The proposer should run \((N, K) = (500, 128), (500, 160),
(1000, 128)\) and extrapolate.

### Comment 4. Inheritance is clean.

v48 inherits from v44 unchanged. No new bugs.

## 5. Feedback to the proposer

For v49:

1. **Run the v48 table** with completed numbers. Report
   \((500, 128), (500, 160), (1000, 128)\).
2. **Compare to the v44 K-sweep prediction**: is \(\Omega(500, 128)\)
   greater than \(\Omega(500, 96) = 1.2814\)? If yes, the \(K\)-sweep
   is still productive; if no, the framework has saturated.
3. **Consider combining v48 + v47** (C–S bounds on \(M[k, k]\)) to
   see if the combination exceeds either alone.
4. **Consolidate process**: stop creating new iterations per
   parameter. Maintain v44 as the rolling "current best framework"
   document.

---

**Final verdict: VALID-but-empty** (v44 run at larger \(K\); no
new math; no new numerical results. Best bound remains v44's
\(C_{6.2}\ge 1.2913\) at \((3000, 96)\)).
