# Verdict — Autocorrelation Inequality 6.2, v6

## 1. Summary of the proof

v6 attempts to extend v5's even-restricted (G) trick to general
admissible \(f\) by introducing **free slack variables** \(s_j \ge 0\)
for \(j = 1,\ldots,2N\) intended to lift the per-index asymmetry
\((p_j - p_{2N+1-j})^{2}\). Two new constraints:

* **(S1)** \(\;s_j \ge Q[j,j] + Q[2N+1-j,\,2N+1-j] - 2 Q[j,\,2N+1-j]\),
* **(S2)** \(\;\sum_j s_j \le 2\,\mathrm{tr}(Q) - 2 \sum_j Q[j,\,2N+1-j]\).

The proof shows that both hold with equality for the rank-one lift
\(Q = pp^{\top}\) and \(s_j = (p_j - p_{2N+1-j})^{2}\), so they are
valid constraints on the relaxation. The proof closes with a numerical
observation that the LP optimum is *still* \(\Omega^{(6)} = 1\)
because the slacks are "phantom" (free, absorbed into the diagonal-\(Q\)
solution), and a v7 plan to switch to the \(F\)-side via Fourier
moments.

## 2. Verdict

**INVALID** (by inheritance, *plus* the new constraints make no
contribution).

The (S1) and (S2) derivations are *correct* — they're valid
constraints satisfied by the rank-one lift of any admissible \(f\).
But:

1. **The new constraints are phantom.** The slacks \(s_j\) are
   introduced as free LP variables with no link to the objective
   \(\Omega\). The diagonal-\(Q\) degeneracy that defeated v1, v2, v3,
   v4 absorbs them effortlessly: with \(p = (1/(2N))\mathbf 1\),
   \(Q = (1/(2N))I\), the constraints (S1) require \(s_j\ge 1/N\),
   (S2) requires \(\sum s_j\le 2\), and these are simultaneously met
   by \(s_j = 1/N\) (giving \(\sum = 2\) exactly). So v6's constraints
   add zero structural information.
2. **The inherited v0 band-constraint hand-wave is still unresolved**,
   for the 6th iteration in a row.

The proof's "Validity" section says

> For every admissible \(f\), the assignment ... is feasible: it
> satisfies (S1) and (S2) with equality and **all v0/v1 constraints**.

so the v6 claim \(\Omega^{(6)}\le C_{6.2}\) inherits the same v0 band
issue.

## 3. Verification of the new content

### (S1). ✓
For the rank-one truth, RHS
\(= p_j^{2} + p_{2N+1-j}^{2} - 2 p_j p_{2N+1-j} = (p_j - p_{2N+1-j})^{2}\),
which is \(\ge 0\). Setting \(s_j := (p_j - p_{2N+1-j})^{2}\) achieves
equality. ✓

### (S2). ✓
For the rank-one truth, summing the equality form of (S1):
\[
\sum_j (p_j - p_{2N+1-j})^{2}
\;=\; \sum_j (p_j^{2} - 2 p_j p_{2N+1-j} + p_{2N+1-j}^{2})
\;=\; 2\!\sum_j p_j^{2} - 2\!\sum_j p_j p_{2N+1-j}.
\]
The first sum is \(\mathrm{tr}(Q)\). The second sum is
\(\sum_j Q[j,\,2N+1-j]\). So
\(\sum_j s_j = 2\,\mathrm{tr}(Q) - 2\sum_j Q[j,\,2N+1-j]\),
which is the equality form of (S2). ✓

(Note: the symmetry-summation step uses
\(\sum_j p_{2N+1-j}^{2} = \sum_j p_j^{2}\), which is true because the
re-indexing \(j\to 2N+1-j\) is a bijection on \(\{1,\ldots,2N\}\).)

### "Validity of the relaxation". ✓ for the new content; ✗ for the inheritance.
The new (S1), (S2) constraints are correctly proved to hold for the
rank-one lift of any admissible \(f\). But the section also asserts
"all v0/v1 constraints" hold "by their respective derivations", which
inherits the v0 band-constraint chain whose middle equality
\(\frac12\sum_{j+k\in\{m,m+1\}} p_j p_k = \int_{B_m}(f*f)\,dt\) is
still unproven. Same as the previous five iterations.

## 4. Why the new constraints are phantom

### The slacks are unanchored to anything that controls \(\Omega\).

v5's (G) constraint
\(Q[j,2N+1-j] = Q[j,j]\) succeeded because it was an *equality* on
the LP variables: it forced the off-diagonal entry of \(Q\) to equal
the diagonal entry, which is *non-trivial* for any \(Q\) with nonzero
trace. The diagonal mixture \(Q = (1/(2N))I\) violated (G) by the
maximum possible amount.

v6's (S1) is an *inequality*: \(s_j \ge \text{(quadratic in } Q)\),
where \(s_j\) is a free LP variable. There is no constraint of the
form "\(s_j \le \text{something involving } \Omega\)" or
"\(\sum_j s_j\) appears anywhere in the objective." So the LP simply
sets \(s_j\) to whatever value is convenient, and (S1) is silently
trivially satisfied. The constraint encodes the *fact* that
asymmetry is nonneg, but doesn't *use* it.

To make the slacks meaningful, the proof would need to:

* Bound \(\sum_j s_j\) **from above** by something *not involving the
  free \(s_j\)* (e.g., by a quantity computable from \(p\) only). This
  would force the asymmetry to be small.
* And/or link \(\sum_j s_j\) to \(\Omega\) via an inequality of the
  form "if \(f\) is asymmetric, then \(\|f*f\|_\infty\) is bounded
  below by something."

Neither is done. (S1) and (S2) on their own provide *no* coupling.

### The proof's diagnosis is correct, repeated for the 6th time.

> The fundamental obstruction is that the moment cone of admissible
> \(f\)'s is the same as the convex hull of single-\(f\) moments,
> which includes mixtures of basis-vector measures. Any convex
> relaxation that uses only \((p, Q, \text{linear or PSD lifts})\) on
> the \(f\)-side cannot strictly improve on \(\Omega = 1\).

This is the same diagnosis as v3's "Numerical observation" about
Lasserre levels, v4's diagonal-\(Q\) analysis, and v5's (G) construction
(which fixed it for *even* \(f\) only). v6 is the 6th iteration to
hit this wall and the 5th failure (only v5 sidestepped it, by
restricting to even \(f\)).

The diagnosis is *correct*. What's missing is acting on it: stop
adding constraints to the \((p, Q)\)-side LP and *actually* switch to
the \(F\)-side or the dual side, as v4 (Section 5 last bullet) and v5
(closing v6 plan) and now v6 (closing v7 plan) have all promised.

## 5. Comments

### Comment 1. v6 doesn't execute v5's stated plan.

v5's closing paragraph said:

> v6 will attempt this distinction by working with both \(f\) and
> \(\check f(x) = f(-x)\) simultaneously and using the fact that
> \(\int f(x) f(-x)\,dx = \frac{1}{2}(\int f^{2} - \int(f - \check f)^{2}/2)\),
> which couples the evenness defect to a quantity that the relaxation
> can constrain.

(I noted in the v5 verdict that the parallelogram identity has a
factor-of-2 error: the correct form is
\(\int f \check f = \int f^{2} - \frac12 \int(f - \check f)^{2}\).)

The *plan* was to use this identity to *couple* the evenness defect
\(\int(f - \check f)^{2}\) to the symmetric inner product
\(\int f \check f = (f * \check f)(0)\), and thereby derive a constraint
that ties the asymmetry to a quantity controlled by \(\Omega\).

**v6 does not execute this plan.** It introduces the slack variables
\(s_j\) but never derives the coupling. There is no constraint of the
form "\(\sum_j s_j \le c \cdot \Omega + d\)" or
"\(\sum_j s_j + \text{(symmetric inner product)} = \int f^{2}\)" that
would let the LP propagate the asymmetry into a bound on \(\Omega\).
The slacks just sit there.

Without the coupling, (S1) and (S2) are vacuous. v6 needed to
finish the v5 plan, and doesn't.

### Comment 2. The v7 plan (\(F\)-side) is the right structural move.

The v7 plan is:

> working entirely on the \(F\)-side with constraints that exclude
> uniform \(F\), which is the approach v7 will explore via Fourier
> moments of \(F\) tied to \(f\) via the squared identity.

This is the same as the 6.3 v4 autocorrelation Toeplitz construction,
and the same as White's Section 4 LP for the Erdős problem. It is the
*right* move and should have been v6 (or even v4). The repeated
"primal-side moment lifts" of v3, v4, v6 have all bounced off the
same wall, and v6's diagnosis correctly identifies it.

**Concrete v7 plan I'd suggest:**

1. Discretize \(F = f * f\) on \([-1/2, 1/2]\) into \(4N\) cells of
   width \(L\); let \(W_m = \int_{B_m} F\,dt/L\) be the cell averages.
   The objective is \(\Omega = \max_m W_m\) (which is \(\le \|F\|_\infty\),
   so a *lower* bound on \(\Omega\) is the right direction for the
   inf-problem).
2. Express \(W_m\) as a quadratic functional of \(p\) (the cell
   averages of \(f\)) via the convolution identity, and lift to
   moments \(y_\alpha\) as in v3's Lasserre setup.
3. **Crucially:** add the *Fourier-side constraint*
   \(\hat F(k) = \hat f(k)^{2}\) (which is exactly the rank-one
   identification on the Fourier side). This is *not* a constraint
   that the basis-vector mixture satisfies, because the basis-vector
   mixture corresponds to *different* \(f\)'s with different
   \(\hat f(k)\)'s, so the average \(\hat F(k)\) is not the *square*
   of an average \(\hat f(k)\).
4. Take a low-degree truncation \(K\) and add the linear constraints
   (Re/Im of \(\hat F(k)\)) to the LP.

This is the natural F-side analog of 6.3 v4's autocorrelation
Toeplitz lift. It has at least a chance of breaking the
\(\Omega = 1\) trap, because the Fourier-side rank-one identification
distinguishes single \(f\)'s from mixtures.

### Comment 3. v6 inherits the v0 band hand-wave for the 6th time.

I've been flagging this since the v1 verdict, and it's still not
addressed. **The right time to fix this is now**, *before* v7's
\(F\)-side switch — because the \(F\)-side LP is going to need its
own band/window constraints, and inheriting v0's hand-wave to a new
framework would propagate the issue further.

The fix sketch (from the v3 and v4 verdict feedback): add cell-supremum
variables \(M_j \ge p_j/L\) (which is the cell-average lower bound on
the supremum) and use them in a *correct* upper bound on
\(\int_{B_m}(f*f)\) via cell-rectangle areas. This makes the band
constraint rigorous and removes the inheritance from every future
verdict.

Three lines of LP-side code; one paragraph of proof. Worth doing.

## 6. Feedback to the proposer

Five iterations of "throw more constraints at the \((p, Q)\)-side LP"
is enough. The diagnosis has been correct since v3, the failure mode
has been identical since v3, and the structural fix has been
identified since v4. **Stop iterating on the primal side. Switch to
the \(F\)-side / Fourier side / dual side, whichever is fastest.**

Concrete recommendations for v7:

1. **Drop v6 entirely.** Phantom slacks add nothing; their absence
   would not change the LP optimum. Keep v5's (G) for the even case
   as a sanity check (and to confirm the LP technology works), but
   don't try to extend (G) to the general case via slacks.
2. **Build the \(F\)-side LP** as sketched in Comment 2. The key
   constraint is the Fourier-side \(\hat F(k) = \hat f(k)^{2}\) (the
   rank-one identification). This is the same idea as 6.3 v4's
   autocorrelation Toeplitz lift, but transposed to the 6.2 problem.
3. **Fix the v0 band-constraint hand-wave** before doing v7. Six
   inheritance verdicts is enough.

If you do all three, v7 should produce \(\Omega^{(7)} > 1\) for the
first time in the project, which would be a real numerical
improvement on the SOTA.

---

**Final verdict: INVALID.**
