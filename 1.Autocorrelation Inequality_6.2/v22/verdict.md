# Verdict — Autocorrelation Inequality 6.2, v22

## 1. State of the file

v22 attempts to tighten the v14 cell-bounds on \(a_k, b_k\) by
replacing the loose worst-case
\(\sum_j p_j(\min\cos)_{j,k}\le a_k\le \sum_j p_j(\max\cos)_{j,k}\)
with a *midpoint quadrature plus error*:
\[
\Bigl|\,a_k - \sum_j p_j\cos(2\pi k x_j)\,\Bigr|\;\le\; \mathrm{err}_k,
\]
with the proposed \(\mathrm{err}_k = \pi^{2} k^{2} L_f^{2}/2\) (the
standard midpoint-rule error for *smooth* integrands).

The proposer **catches their own mistake** before submitting: they
verify that this error bound is too tight for non-smooth (e.g.
delta-like) admissible \(f\), and explicitly mark v22 as INVALID
with a v23 plan.

## 2. Verdict

**INVALID** (self-INVALID, correctly diagnosed).

The proposer is right. The midpoint rule's \(O(L^{2})\) error
formula assumes the integrand has bounded second derivative — it
does not apply to the LP integrand \(f \cdot \cos(2\pi k x)\) when
\(f\) can concentrate (e.g., delta-like). For non-smooth \(f\), the
error is \(O(L_f)\), not \(O(L_f^{2})\), and v22's constraint is
too tight by a factor of \(\sim 1/L_f\).

This is the same kind of "tried a constraint, tested on a known
admissible \(f\), found a counterexample, marked self-INVALID"
workflow as 6.2 v8. **The right kind of failure.**

## 3. Verification of the self-INVALID

I checked the proposer's counterexample.

### The delta example.

For \(f = \delta_{x^{*}}\) (concentrated at \(x^{*}\in I_j\)) with
\(p_j = 1\):

* **Exact**: \(\int_{I_j}f(x)\cos(2\pi k x)\,dx = \cos(2\pi k x^{*})\).
* **Midpoint approximation**: \(p_j\cos(2\pi k x_j) = \cos(2\pi k x_j)\),
  where \(x_j\) is the cell midpoint.

The error:
\[
\big|\cos(2\pi k x^{*}) - \cos(2\pi k x_j)\big|.
\]

By the mean-value theorem (or directly bounding the derivative
\(|\cos'| \le 1\)):
\[
|\cos(a) - \cos(b)| \;\le\; |a - b|.
\]
With \(a = 2\pi k x^{*}\), \(b = 2\pi k x_j\), and
\(|x^{*} - x_j|\le L_f/2\) (since \(x^{*}\in I_j\)):
\[
|\cos(2\pi k x^{*}) - \cos(2\pi k x_j)| \;\le\; 2\pi k\cdot L_f/2 \;=\; \pi k L_f.
\]

So the worst-case error is \(\pi k L_f\), which is \(O(L_f)\), *not*
\(O(L_f^{2})\). The proposer's bound \(\pi^{2} k^{2} L_f^{2}/2\) is
too tight by a factor of \(\sim 1/(\pi k L_f) = O(L_f^{-1})\) at
fixed \(k\).

### Numerical sanity check.

For \(N = 8\), \(L_f = 1/32\), \(k = 1\):
* True worst-case error: \(\pi\cdot 1\cdot 1/32 = \pi/32 \approx 0.0982\).
* v22's allowed error: \(\pi^{2}\cdot 1\cdot (1/32)^{2}/2 = \pi^{2}/2048 \approx 0.00482\).
* Ratio: \(0.0982/0.00482 \approx 20\). ✓

So v22's constraint is too tight by \(\sim 20\times\) at this
parameter setting, and excludes the delta-concentrated \(f\) (or
its \(L^{1}\) approximations).

### Why the standard midpoint error doesn't apply.

The standard midpoint rule error
\(|f(x_j) \cdot L_f - \int_{I_j}f(x)\,dx|\le L_f^{3}\sup|f''|/24\)
*requires* \(f\) to have a bounded second derivative. In v14's LP,
\(f\) is only \(L^{1}\) (a nonneg density), and can concentrate to
delta-like behaviour where \(f''\) does not exist (in the classical
sense) or is unbounded. The midpoint rule does not give an
\(O(L^{2})\) bound for non-smooth integrands.

The *correct* worst-case error for non-smooth \(f\) is \(O(L)\),
which is exactly the loose v14 bound \(\sum p_j(\min\cos)_{j,k}\le a_k\le \sum p_j(\max\cos)_{j,k}\)
in disguise — the v14 bound's slack is also \(O(L)\) per cell.

So replacing v14's worst-case bound with the proposer's "midpoint +
\(O(L^{2})\) error" doesn't actually tighten anything *and* breaks
rigor. v14's bound is the right thing.

### The numerical observation.

The proposer reports \(\Omega^{(22)} = 1.124\) at \((N, K) = (32, 32)\),
which is *higher* than v14's \(\sim 1.088\) at the same
parameters. This confirms that the v22 constraint is genuinely
*tighter* (it gives a larger LP optimum, which is closer to the
true \(C_{6.2}\)). But because the constraint is *non-rigorous*
(too tight), the larger bound is *not certifiably* a lower bound
on \(C_{6.2}\).

It's possible that \(\Omega^{(22)} = 1.124 > C_{6.2}\), in which
case v22 is wrong. The proposer is right to mark it INVALID and
not claim the bound.

## 4. Comments

### Comment 1. The self-checking is excellent.

This is the second time in 6.2 (after v8) that the proposer has
caught their own mistake before submitting. The pattern:

1. Try a constraint.
2. Test on a specific admissible \(f\) (delta-concentrated, in
   this case).
3. Compute the actual error.
4. Compare to the proposed bound.
5. If the proposed bound is violated, mark INVALID.

This is the *right* workflow, and the verdict instructions
explicitly ask for it. Each self-INVALID is more valuable than ten
"VALID, no improvement" iterations because it teaches the project
exactly what doesn't work and why.

### Comment 2. The structural lesson.

The lesson from v22: **the LP integrand smoothness is the binding
constraint on the discretization tightness.** As long as the LP
allows \(f\) to be unboundedly concentrated (which is true under
the standard \(f\ge 0, \int f = 1, \operatorname{supp}(f)\subset[-1/4,1/4]\)
assumptions), the discretization error is \(O(L)\), not \(O(L^{2})\).

To get a tighter discretization, you'd need to *restrict* to \(f\)
with bounded \(L^{\infty}\) norm (\(\|f\|_{\infty}\le M\) for some
\(M\)), which would give a tighter error in terms of \(M\). But this
restricts the bound to \(C_{6.2}^{M}\), the bounded-\(M\) version of
\(C_{6.2}\), which is *strictly larger* than \(C_{6.2}\) (since
\(C_{6.2}\) takes the inf over all \(f\)).

The proposer's v23 plan correctly identifies this trade-off:

> Adding a constraint that \(f\) is bounded ... would justify a
> higher-order quadrature error. However, \(\|f\|_{\infty}\) is *not
> bounded above* on the admissible set ... Option 2 doesn't work
> for arbitrary admissible \(f\), but it could work if we *restrict*
> to \(f\) bounded by a known \(M\) — at the cost of bounding only
> \(C_{6.2}^{M}\), which \(\le C_{6.2}\) in general.

This is the right framing. v23's restricted bound is *not* a bound
on the original \(C_{6.2}\) — it's a bound on the related
constant \(C_{6.2}^{M}\). For a *true* improvement on the SOTA,
the proposer would need to argue that \(C_{6.2}^{M}\) approximates
\(C_{6.2}\) for some computable \(M\), which is a separate analytical
question.

### Comment 3. The numerical bounds are interesting even though invalid.

The v22 bound 1.124 at \((32, 32)\) — even though *not certified*
as a lower bound on \(C_{6.2}\) — is informative:

* It suggests that the v14 framework can *in principle* reach
  bounds well above the v14 plateau if the cell-bound on \(a_k\)
  is tightened.
* It suggests that the *non-rigorous* (heuristic) version of the
  LP gives bounds approaching the SOTA 1.28.
* It tells the proposer that the discretization slack is the
  binding bottleneck, not the SDP framework itself.

So v22's failure is *informative* even though it's not a valid
result.

### Comment 4. The right next move (v23).

The v23 plan is to add a constraint \(\|f\|_{\infty}\le M\) and
bound \(C_{6.2}^{M}\). This is rigorous *for the restricted
problem*, but doesn't directly improve the bound on \(C_{6.2}\).

To make v23 useful for the original problem, the proposer should
also:
1. Show that \(C_{6.2}^{M} \ge \) some explicit lower bound on
   \(C_{6.2}\) for sufficiently large \(M\), or
2. Take the inf over \(M\) (which recovers \(C_{6.2}\) at \(M\to\infty\),
   but the bound at finite \(M\) is only valid for \(C_{6.2}^{M}\)).

This is non-trivial. An alternative: prove that the *infimum-attaining
\(f\)* in \(C_{6.2}\) has bounded \(\|f\|_{\infty}\), so the
restriction is WLOG. This is plausible (the infimum is probably
achieved by some "smooth" \(f\)) but requires proof.

## 5. Feedback to the proposer

Excellent self-checking. The verdict on v22 is INVALID, but the
*workflow* is exactly right. Each self-INVALID iteration teaches
the project exactly what doesn't work.

For v23:
1. **Implement the bounded-\(M\) version** as planned. The bound
   \(C_{6.2}^{M}\) is a *different* constant from \(C_{6.2}\), and
   the proposer should be careful to label it correctly.
2. **Investigate whether the infimum in \(C_{6.2}\) is achieved by
   bounded \(f\).** If yes, then \(C_{6.2}^{M} = C_{6.2}\) for
   large enough \(M\), and v23's bound is also a bound on the
   original. If no, the gap between \(C_{6.2}^{M}\) and \(C_{6.2}\)
   needs to be quantified.
3. **As an alternative**, try a *moment-lifted* version of the
   discretization: instead of bounding \(\|f\|_{\infty}\), bound
   \(\int f^{2}\) (which is an LP-side moment) by an upper-bound
   parameter. This is closer in spirit to the v14 framework and
   doesn't require restricting the admissible class.

Either way, v23 should produce a *legitimately rigorous* bound,
not a heuristic one.

---

**Final verdict: INVALID** (self-INVALID, correct diagnosis).
