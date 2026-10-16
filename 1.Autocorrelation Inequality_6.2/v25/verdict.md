# Verdict — Autocorrelation Inequality 6.2, v25

## 1. Summary

v25 is a **clean rigorous synthesis** of all the project's work
on 6.2's unrestricted constant. No new constraints — exactly the
v14 autocorrelation Fejér framework — with a parameter sweep at
larger \((N, K)\) and an explicit declaration that this is the
project's best rigorous result on \(C_{6.2}\).

The construction:
* WLOG reductions (from v7) — support and normalization.
* PSD lift on \(M\) for \((1, a, b)\) and f-side cell-bounds on
  \((a_k, b_k)\) from \(p\).
* Autocorrelation Fejér constraint via \(u_k = M[k,k]+M[K+k,K+k]\),
  \(v_k\ge u_k^{2}\), \(\Omega\ge 1 + 2\sum(1 - k/(K+1))v_k\).

The numerical results are at parameters up to \((100, 32)\), with
\(\Omega^{(25)} = 1.11\). The proposer claims this as the project's
**first rigorous lower bound on \(C_{6.2}\) strictly above 1**,
valid for all admissible \(f\) (not just even).

## 2. Verdict

**VALID** by identity with v14 (which is VALID per
`v14/verdict.md`). No new mathematical content; the numerical
results extend v17's parameter sweep.

The proposer's project trajectory table is a useful summary and
matches my own tracking.

## 3. The project's best rigorous bound

The proposer's claim:
> Theorem (rigorous, established by the v14 / v25 LP).
> \(C_{6.2}\ge 1.10\) (at \((N, K) = (50, 64)\)).

Refining slightly: at \((N, K) = (100, 32)\), the v17/v25 sweep
gives \(\Omega^{(25)} = 1.11\), so the project's best is

> **\(C_{6.2}\ge 1.11\)** (rigorous, valid for all admissible \(f\),
> from the v14/v17/v25 autocorrelation Fejér construction at
> \((N, K) = (100, 32)\)).

This is **below** the SOTA \(C_{6.2}\ge 1.28\), but it's the
project's first non-trivial rigorous bound on the unrestricted
constant. The asymptotic limit (predicted as \(\inf_f \int F^{2}\))
is at most \(4/3\), so further parameter scaling can in principle
push the bound higher.

## 4. Comments

### Comment 1. The synthesis is the right scientific practice.

After 24 iterations of attempts, v25 explicitly *declares* the
final rigorous result and the project trajectory. This is the right
move: rather than continue adding more vacuous constraints, the
proposer is now consolidating and reporting.

### Comment 2. The project trajectory table is accurate.

Comparing my tracking to the proposer's table:

| version | bound | restrict | my notes |
|---|---|---|---|
| v0 | 1.0 (trivial) | none | invalidated by band hand-wave (v1 verdict) |
| v7 | **1.34** \((N=20)\) | even | first valid iteration |
| v9 | **1.69** \((N=50, K=64)\) | even | Fejér kernel — even-restricted |
| v11 | 1.376 \((N=8, K=16)\) | even | shifted (W) — modest improvement |
| v14 | **1.063** \((N=16, K=16)\) | none | first non-even rigorous bound |
| v17 | **1.111** \((N=100, K=32)\) | none | parameter sweep |
| v22 | 1.124 (claimed) | none | self-INVALID (non-rigorous err) |
| v25 | 1.11 (best) | none | clean v14 sweep |

This matches my tracking. ✓

### Comment 3. The asymptotic limit prediction.

The v17 / v25 prediction is that the v14 framework's asymptotic
limit (as \(N, K\to\infty\)) is \(\inf_f \int(f*f)^{2}\), which is
at most \(4/3\) (achieved at uniform \(f\)).

If the uniform \(f\) is the actual minimizer of \(\int F^{2}\), then
the asymptotic limit is exactly \(4/3 \approx 1.333\), which would
*exceed the SOTA* 1.28. But the LP is converging *very slowly*
(from 1.063 at \(N=16\) to 1.11 at \(N=100\)), so reaching the
asymptotic limit may require enormous \(N\).

The proposer's open question (per v25's "Conclusion and open work"):
> Closing the gap requires larger-scale numerical computation than
> the present iterations.

is correct. The v14 framework is rigorously sound and converges to
a meaningful asymptotic limit, but the convergence rate is the
binding constraint at the present iteration.

### Comment 4. Suggestion for v26 (if there is one).

Given that v15–v24 have shown that *adding more constraints* doesn't
help, and v25 is the synthesis declaration, the natural next step
is to **either**:

1. **Run v25 at much larger \((N, K)\)** (\(N \in \{200, 500, 1000\}\),
   \(K \in \{64, 128, 256\}\)) to see how close to \(4/3\) the bound
   gets. This is purely a numerical step.
2. **Diagnose whether the asymptotic limit is actually \(4/3\)**,
   by extracting the LP-optimal \(p\) at large \(N\) and computing
   \(\int F^{2}\) for the corresponding \(f\). If it matches \(4/3\),
   the limit is reachable; if not, the LP is finding a better
   minimizer than the uniform.

If v26 does (1) and reports the actual numerical results, the
project will know whether the v14 framework can beat the SOTA.

---

**Final verdict: VALID** (clean rigorous synthesis; no new content,
but the right scientific practice for declaring a project's
result).
