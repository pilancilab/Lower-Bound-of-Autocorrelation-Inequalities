# Verdict — Autocorrelation Inequality 6.2, v49

## 1. Summary

v49 is another **parameter-push iteration**, this time for v45
(= v44 + Bochner Toeplitz at degree \(K\)) at larger \(K = 128, 160\)
at moderate \(N\). No new mathematical content. Running in parallel
with v48 (v44 at larger \(K\)).

**No numerical results yet.** The best rigorous bound remains
v44's \(C_{6.2}\ge 1.2913\) at \((3000, 96)\).

## 2. Verdict

**VALID-but-empty** (v45 at larger parameters; no new math; no
new results).

## 3. Comments

### Comment 1. Same process issue as v48.

v48 and v49 are both "run existing framework at larger \(K\)" with
no new math. They should be consolidated into a single "production
sweep" document, not split into per-\(K\) iterations.

### Comment 2. The experiment is sensible.

Comparing v48 (v44 alone) vs v49 (v44 + Bochner) at \(K = 128, 160\)
is a reasonable experiment to see whether the Bochner constraint
becomes more or less binding at higher \(K\). Given v45's
production-sweep result showed Bochner adds only \(\sim 10^{-5}\)
at \(K = 64, 96\), my prediction is that it remains \(\sim 10^{-5}\)
at \(K = 128, 160\) — the truncated Bochner is essentially captured
by the v44 moment lift at any \(K\).

If the Bochner constraint slows the solve prohibitively (as the
proposer worries), then v45 at high \(K\) is strictly dominated by
v44 (tiny gain, large solve-time cost).

### Comment 3. Inheritance is clean.

v49 inherits from v45. No new bugs.

## 4. Feedback to the proposer

1. **Complete the v48 and v49 parameter runs** and tabulate the
   results in a single document.
2. **If v45's Bochner adds \(< 10^{-4}\) at \(K = 128\)**, stop
   running v45/v49 and focus on v44/v48.
3. **Consider v47 (C–S on \(M[k,k]\)) at higher \(K\)** — that
   construction had stronger theoretical motivation for being
   binding.
4. **Update the project summary** with v44's \(1.2913\) as the
   current best and clearly mark v48/v49 as parameter runs.

---

**Final verdict: VALID-but-empty** (v45 at larger parameters;
no new math or results; best bound remains v44's \(1.2913\) at
\((3000, 96)\)).
