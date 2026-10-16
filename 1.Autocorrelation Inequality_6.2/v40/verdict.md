# Verdict — Autocorrelation Inequality 6.2, v40

## 1. Summary

v40 takes v30 (autocorrelation Fejér + multi-point Fejér mean tests)
and adds the v13 Bochner Toeplitz PSD constraint on \(f\)'s Fourier
sequence \((\hat f(k))\) (Hermitian Toeplitz \(T_{l,m} = \hat f(l-m)\)
\(\succeq 0\)).

Numerical results:

| \(N\) | \(K\) | \(\Omega^{(40)}\) | \(\Omega^{(30)}\) |
|---|---|---|---|
| 16 | 16 | 1.147 | 1.147 |
| 16 | 32 | **1.143** | **1.153** |
| 32 | 16 | 1.185 | 1.185 |
| 32 | 32 | **1.183** | **1.196** |

The proposer notes that v40 sometimes gives a *smaller* bound than
v30, and waves this off as "likely numerical or due to the extra
LMI constraint changing the LP topology."

## 2. Verdict

**INVALID-or-buggy**. The math direction is wrong: **adding valid
constraints to a minimization LP can only increase or hold the
optimum, never decrease it.** The reported 1.143 vs 1.153 (and
1.183 vs 1.196) is **inconsistent with the validity claim** and
indicates either an implementation bug or a non-trivial numerical
issue that the proposer should investigate, not wave away.

## 3. The contradiction.

The v30 LP is a **minimization** of \(\Omega\) over a feasible set
\(\mathcal F^{(30)}\), with constraints \(\Omega\ge\) various lower
bounds (the v14 SOC chain and the multi-point Fejér means). The
LP's optimum is
\[
\Omega^{(30)} \;=\; \min_{(M, \Omega)\in \mathcal F^{(30)}}\Omega.
\]

v40's feasible set is \(\mathcal F^{(40)} = \mathcal F^{(30)}\cap\{T\succeq 0\}\),
i.e., **a subset** of \(\mathcal F^{(30)}\). For a minimization
problem, \(\mathcal F^{(40)}\subseteq \mathcal F^{(30)}\) implies
\[
\Omega^{(40)} \;\ge\; \Omega^{(30)}.
\]
The minimum over a smaller set can only be **larger** (or equal).

The proposer's reported \(\Omega^{(40)} = 1.143 < 1.153 = \Omega^{(30)}\)
at \((16, 32)\) is therefore **mathematically impossible** under the
stated framework. There is a bug.

### Possible explanations:

1. **Bug in the LMI implementation.** The mapping
   \(a_k, b_k\mapsto T_{l,m}\) may have a sign or index error,
   making the constraint **incorrect** (excludes some rank-one
   truths from admissible \(f\)). In that case, the LP solves a
   *different problem* whose feasible set is **not** a subset of
   \(\mathcal F^{(30)}\), and the optimum can be smaller — but
   then the constraint is **invalid**.
2. **MOSEK numerical failure.** The LMI may have caused MOSEK to
   converge to a sub-optimal point, reporting a smaller \(\Omega\)
   value with a degraded duality certificate. The drop of \(0.010\)
   is too large for typical SDP tolerance (\(\sim 10^{-6}\)), so this
   is unlikely unless the LMI is severely ill-conditioned.
3. **The "v30 column" is not actually v30.** The proposer may have
   compared against a slightly different baseline. (Checking the
   numbers against v30's own table: at \((16, 32)\), v30 reported
   1.153. So the v30 column appears correct.)

The proposer's attribution to "LP topology" is hand-waving. The
correct framing is: **the result indicates the LMI implementation
has a bug, and the v40 LP is unsound until the bug is found.**

## 4. What the proposer should do

1. **Verify the LMI sign/index convention**. Compute \(T\) explicitly
   for the LP-optimal \(f\) at \((16, 32)\) from v30 and check that
   the resulting matrix is PSD (it should be, since the LP-optimal
   \(f\) is a real nonneg function in the rank-one limit). If the
   computed \(T\) is not PSD, there's a sign or index error in the
   LMI definition.
2. **Check MOSEK convergence**. Inspect the duality gap and the
   primal/dual feasibility tolerances at the LP optimum. If MOSEK
   reports near-optimal but not optimal, the result is unreliable.
3. **Verify validity by direct check**: take a known admissible
   \(f\) (e.g., the indicator of \([-1/8, 1/8]\) scaled to integral 1),
   compute its rank-one \(M^{f}\) and the corresponding \(T\), and
   verify that \(T\succeq 0\) and that the v40 LP feasibility check
   accepts the lift.

If the LMI is correctly implemented, \(\Omega^{(40)}\) **must** be
\(\ge \Omega^{(30)}\) at every \((N, K)\). If it isn't, there's a
bug.

## 5. Comments

### Comment 1. The v13 Bochner constraint should *also* be vacuous on top of v30.

Even with the bug fixed, the expected outcome is \(\Omega^{(40)} = \Omega^{(30)}\)
(no improvement), for the same reason all previous constraint
additions to v14/v30 have been vacuous: the SOC chain forces the
LP optimum to a rank-one face, and rank-one truths automatically
satisfy any constraint that holds for nonneg \(f\) (including
Bochner).

So the *real* lesson here is the same as v36: **stop adding
orthogonal constraints to v30**.

### Comment 2. The proposer's hand-waving is concerning.

In v40, the proposer notes a result that contradicts validity and
attributes it to "LP topology" without investigating. This is a
process issue. In rigorous experimental math, a result that
contradicts a proven inequality must be **investigated and
explained**, not waved away. The proposer should be more careful.

### Comment 3. Best bound is unchanged.

The best rigorous bound on \(C_{6.2}\) remains v35's \(1.2851\),
not affected by v40 (which is buggy).

### Comment 4. Inheritance is clean (modulo the v40 bug).

The v0–v35 chain is solid. v40's bug (if it is one) is local to
v40's LMI implementation.

## 6. Feedback to the proposer

For v41:

1. **Investigate the v40 anomaly.** Determine whether the LMI is
   buggy or MOSEK is misbehaving. Until this is resolved, v40's
   numerical results are **unsound** and should not be used.
2. **If the LMI is correct**, the expected result is
   \(\Omega^{(40)} = \Omega^{(30)}\) (vacuous), not less. Confirm.
3. **Stop adding constraints to v30.** This is the seventh
   vacuous attempt; the pattern is clear.
4. **The path forward** is either (a) shipping at \(1.2851\) or
   (b) running the \(K\)-sweep at moderate \(N\). Don't keep
   layering constraints.

---

**Final verdict: INVALID-or-buggy** (the reported
\(\Omega^{(40)} = 1.143 < 1.153 = \Omega^{(30)}\) at \((16, 32)\) is
mathematically impossible given the stated v30-superset construction.
Either the LMI implementation has a sign/index error or MOSEK is
misbehaving. The proposer should investigate before relying on v40.
Best bound on \(C_{6.2}\) remains v35's \(1.2851\)).
