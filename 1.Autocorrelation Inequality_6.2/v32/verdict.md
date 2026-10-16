# Verdict — Autocorrelation Inequality 6.2, v32

## 1. Summary

v32 takes v30 (autocorrelation Fejér + pointwise Fejér mean tests on
\(F\)) and adds the v11/v7 shifted windowed (W) integral
inequalities on the \((p, Q)\) lift. The combination is rigorously
valid (each constraint is individually valid for any rank-one truth)
but produces **no numerical improvement** at \((N, K) = (16, 16)\):

| Setting | Bound |
|---|---|
| v30 at \((16, 16)\) | 1.147 |
| **v32 at \((16, 16)\)** | **1.147** (identical) |

The (W) constraints are **not binding** at the v30 LP optimum.

## 2. Verdict

**VALID** (rigorous, but vacuous tightening).

## 3. Verification

Each ingredient was previously verified:
* v14 autocorrelation Fejér chain — verified in v14's verdict.
* v30 pointwise Fejér mean tests — verified in v30's verdict.
* v7/v11 shifted (W) windowed inequalities — verified in v7's
  verdict (these are integral inequalities of the form
  \(\int_{W} F\le c_{W}\) for various window sets \(W\)).
* The \((p, Q)\) PSD lift, mass, marginals, \(Q\ge 0\) — verified in
  the v0–v7 chain.

The conjunction is valid: the rank-one truth from any admissible
\(f\) satisfies *every* one of these constraints simultaneously, so
the LP feasible set still contains the rank-one truth, and
\(\Omega^{(32)}\le C_{6.2}\). ✓

The numerical observation \(\Omega^{(32)} = 1.147 = \Omega^{(30)}\)
at \((16, 16)\) shows that the (W) constraints lie on the *interior*
of the v30 feasible region (not binding at the optimum).

## 4. Comments

### Comment 1. The proposer did not take the v31 verdict's primary recommendation.

The v31 verdict closed by recommending a **\(K\)-sweep** at fixed
\(N\) (e.g., \((50, 100), (50, 200), (100, 100), (100, 200)\)), to
test v29's Nyquist heuristic. Instead, v32 added an orthogonal
constraint family (the (W) windows), which the v28-class history
strongly suggested would not help (v15–v28 showed that "extra"
constraints layered on top of v14 are essentially always vacuous,
because v14's SOC chain already forces the LP to behave as if
rank-one at the optimum).

The proposer should run the \(K\)-sweep before adding more
constraint types.

### Comment 2. Why (W) is dominated.

The (W) constraints are integral inequalities
\(\int_{W} (f * f)(t)\,dt\le c_{W}\) for various window sets \(W\).
At the rank-one truth, these are *averages* of \(F\) over windows
of width comparable to \(L_f\). The pointwise Fejér mean tests are
*averages* of \(F\) over windows of width \(\sim 1/(K+1)\) (the
Fejér kernel width). For \(K = 16\), \(1/(K+1) \approx 0.06\),
which is comparable to or smaller than typical (W) window widths.

So the v30 pointwise Fejér tests are at least as tight as the (W)
windowed averages, and the (W) constraints add no new information
at the v30 LP optimum.

This is consistent with the proposer's observation that
\(\Omega^{(32)} = \Omega^{(30)}\).

### Comment 3. Inheritance is clean.

v32 inherits from v30 + v7/v11. No new bugs.

## 5. Feedback to the proposer

v32 is a rigorous but vacuous iteration. The path forward is the
**\(K\)-sweep** at fixed \(N\), per the v31 verdict's recommendation:

1. \((N, K) = (50, 50), (50, 100), (50, 200)\) — test the Nyquist
   prediction at moderate \(N\).
2. \((N, K) = (100, 100), (100, 200), (100, 400)\) — same at larger
   \(N\).
3. **Goal**: see whether \(K \sim 4N\) gives a substantially better
   bound (e.g., \(\ge 1.28\)) than \(K = 32\).

If the bound exceeds 1.28 at any setting, the project matches the
SOTA. If it plateaus below 1.28, the v30 construction has a true
ceiling, and the proposer should look for a fundamentally new
constraint type (e.g., higher-order SOS hierarchies on \(f\)).

**Stop adding orthogonal constraint families** until the parameter
sweep is complete. v15–v28's history shows that adding constraint
types is the wrong way forward — every one of them was vacuous on
top of v14.

---

**Final verdict: VALID** (rigorous, but the (W) constraints are
not binding at the v30 LP optimum and produce no improvement.
\(\Omega^{(32)} = \Omega^{(30)} = 1.147\) at \((16, 16)\). Best
bound on \(C_{6.2}\) remains v31's \(1.267\) at \((2000, 32)\)).
