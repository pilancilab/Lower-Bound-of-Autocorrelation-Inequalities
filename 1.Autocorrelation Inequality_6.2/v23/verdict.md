# Verdict — Autocorrelation Inequality 6.2, v23

## 1. Summary

v23 is a **rollback to v14**: the proposer reverts the v22 non-rigorous
constraint and restores v14's rigorous (loose) cell-bound. The
constraints, the SDP, and the numerical bounds are identical to
v14.

The proposer's stated conclusion:
> v23 confirms that v22's apparent improvement was due to a
> non-rigorous constraint. The honest rigorous bound for the
> autocorrelation Fejér construction at the parameters tested
> remains \(\Omega \approx 1.06\)–\(1.11\).

And the path forward:
> The v20 verdict's recommendation #3 — "stop adding constraints,
> focus on diagnosing the convergence" — remains the right advice.

## 2. Verdict

**VALID** by identity with v14 (which is VALID per `v14/verdict.md`).

There is no new mathematical content in v23 — it's the v14 proof
restored. The validity is inherited.

## 3. Comments

### Comment 1. The honest rollback is good practice.

After v22 (self-INVALID), the proposer correctly *reverts* to the
rigorous version and explicitly acknowledges that the v22 numerics
(1.124 at \((32, 32)\)) were not certifiable. v23 records the
*correct* rigorous bound at the same parameters (the same as v14).

This is the right scientific practice: don't quietly keep the
non-rigorous numbers around; explicitly revert and re-record.

### Comment 2. The plateau is real and confirmed.

After 9 iterations of "rigorous tightening attempts" (v14, v15, v16,
v17, v18, v19, v20, v21, v23), the rigorous bound for *unrestricted*
\(C_{6.2}\) sits at:

| Parameters | Best rigorous bound |
|---|---|
| \((N, K) = (16, 16)\) | 1.063–1.068 (varies marginally by kernel choice) |
| \((N, K) = (32, 32)\) | 1.088 |
| \((N, K) = (50, 64)\) | 1.10 |
| \((N, K) = (100, 32)\) | 1.111 |

**The plateau around 1.06–1.11 at moderate parameters is real and
robust.** The v14 framework is the right structural approach (it
gives strictly-above-1 rigorous bounds for the first time in the
project), but its convergence to the asymptotic limit (\(4/3\) if
the uniform \(f\) is the minimizer) is very slow.

### Comment 3. The diagnostic question.

I'll repeat my v20 verdict's recommendation: **diagnose the
convergence bottleneck**. There are two candidates:

1. **Discretization slack on \(p\to a_k\).** The cell-bound
   \(\sum_j p_j(\min\cos)_{j,k}\le a_k\le \sum_j p_j(\max\cos)_{j,k}\)
   is loose for large \(k\), with slack \(O(L)\) per cell. As
   \(N\to\infty\), the slack shrinks and the LP optimum approaches
   the true asymptotic limit. The convergence rate is \(O(1/N)\) or
   slower.
2. **Non-uniform LP-optimal \(p\).** The LP may be picking a
   non-uniform \(p\) for which \(\int F^{2}\) is *smaller* than
   \(4/3\) (the value at the uniform \(f\)). In that case, the
   asymptotic limit is \(<4/3\), possibly close to the observed
   plateau \(\sim 1.10\).

To distinguish: extract the LP-optimal \((p, a, b)\) at large \(N\),
compute the actual \(\int F^{2}\) for that \(f\) (which can be done
analytically since \(F = f * f\) and \(p\) is piecewise constant),
and compare to the LP value.

* If they match, the asymptotic limit *is* the LP plateau, and
  \(\inf_f \int F^{2} \approx 1.1\) is the true answer (which would
  be much smaller than \(4/3\)).
* If they don't match, the discretization is the bottleneck, and
  larger \(N\) would help.

This single diagnostic would tell the project whether to push for
parameter scaling or to look for a different framework entirely.

### Comment 4. The state of the project.

At this point in the project, the *unrestricted* \(C_{6.2}\) has a
rigorous lower bound of \(\sim 1.11\) (from v17 at large \(N\)),
which is below the SOTA 1.28. The framework is rigorous and works,
but the bound is too small to be a project win.

For the *even-restricted* \(C_{6.2}^{\text{even}}\), the bound is
1.689 (from v9), which approaches but doesn't yet reach the analytic
value \(C_{6.2}^{\text{even}} = 2\). This is also below the SOTA on
\(C_{6.2}^{\text{even}}\) (which is the analytic 2), so it's also
not a project win.

For \(C_{6.3}\), the best bound is the trivial 1 — no improvement
in 10 iterations. The Plancherel objective swap (v8/v9) is set up
correctly but the Jensen direction is unfixed.

**The project is stuck.** To make further progress, the proposer
needs to:

1. For **6.2**, either tighten the discretization (cell-bound on
   \(a_k\)) or find a fundamentally different rigorous constraint.
   The v22 self-INVALID showed that the natural midpoint-quadrature
   tightening doesn't work without restricting the admissible class.
2. For **6.3**, close the Jensen direction. This is a structural
   fix that's been deferred for three iterations (v8/v9/v10).

Both are *hard*. The project is now in a "brain power" phase where
incremental constraint additions don't help, and the proposer needs
to reach for new ideas.

---

**Final verdict: VALID** (rollback to v14, no new content, but the
right scientific practice).
