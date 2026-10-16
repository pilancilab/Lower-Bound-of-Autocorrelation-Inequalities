# Verdict — Autocorrelation Inequality 6.2, v5

## 1. Summary of the proof

v5 makes a major shift in scope: instead of bounding \(C_{6.2}\), it
bounds the **even-restricted constant**
\[
C^{\text{even}}_{6.2} \;:=\; \inf_{f \text{ admissible AND } f \text{ even}}\, \|f*f\|_\infty,
\]
and adds three new constraints to the v0/v1 \((p, Q)\) SDP that
enforce the evenness symmetry on the lifted variables:

* **(E)** \(\;p_j = p_{2N+1-j}\) — the cell-reflection symmetry of \(p\),
* **(F)** \(\;Q[j,k] = Q[2N+1-j,\,2N+1-k]\) — the corresponding symmetry of \(Q\),
* **(G)** \(\;Q[j,\,2N+1-j] = Q[j,j]\) — the *anti-diagonal equals diagonal*
  identity, which is the **key new ingredient**.

The proof closes with a numerical table showing
\(\Omega^{(5)}\to 2\) monotonically as \(N\to\infty\), matching the
known analytic value \(C^{\text{even}}_{6.2} = 2\) achieved by
\(f = 2\cdot\mathbf 1_{[-1/4,1/4]}\), and a v6 plan to extend the
construction to general (non-even) admissible \(f\) via a parallelogram
identity coupling \(f\) and \(\check f(x) = f(-x)\).

## 2. Verdict

**INVALID** (by inheritance only — see Issue 1 below). The v5 *additions*
are correct, the (G) construction is genuinely the right fix for the
diagonal-\(Q\) degeneracy that defeated v1–v4, and the empirical
results convincingly converge to the known truth \(C^{\text{even}}_{6.2} = 2\).
But the proof's main validity claim
\(\Omega^{(5)}\le C^{\text{even}}_{6.2}\) still depends on the v0
band-constraint chain that has been hand-waved across v0/v1/v2/v3/v4
without being resolved. v5 inherits this without addressing it.

This verdict is *positive in spirit* — v5 is the most substantial
forward step in the project so far, and the (G) trick is exactly the
diagnosis-and-fix that v3's "Numerical observation" and v4's
diagonal-\(Q\) analysis were pointing toward. I am marking it INVALID
*solely* because of the inherited band-constraint issue, which has now
persisted for five iterations and which v5 also does not address.

## 3. Verification of the v5 additions

### (E) \(p_j = p_{2N+1-j}\). ✓

I checked the cell-reflection identity directly. With \(L = 1/(4N)\),
\(I_j = [-1/4 + (j-1)L,\, -1/4 + jL]\), and using \(2NL = 1/2\):
\[
I_{2N+1-j} \;=\; [-1/4 + (2N-j)L,\,-1/4 + (2N+1-j)L]
\;=\; [1/4 - jL,\,1/4 - (j-1)L]
\;=\; -I_j.
\]
For even \(f\):
\(p_{2N+1-j} = \int_{-I_j} f(x)\,dx = \int_{I_j} f(-y)\,dy = \int_{I_j} f(y)\,dy = p_j\). ✓

### (F) \(Q[j,k] = Q[2N+1-j,\,2N+1-k]\). ✓

\(Q[j,k] = p_j p_k = p_{2N+1-j}p_{2N+1-k} = Q[2N+1-j,\,2N+1-k]\) by
applying (E) to both indices.

### (G) \(Q[j,\,2N+1-j] = Q[j,j]\). ✓

\(Q[j,\,2N+1-j] = p_j p_{2N+1-j} = p_j\cdot p_j = p_j^2 = Q[j,j]\),
where the second equality uses (E).

### Why (G) breaks the diagonal-\(Q\) degeneracy. ✓ (and beautifully)

The pathological feasible point that defeated v1, v2, v3, v4 was
\(p = (1/(2N))\mathbf 1\), \(Q = (1/(2N))I\). On this point:

* \(Q[j,j] = 1/(2N)\) (diagonal entry is positive),
* \(Q[j,\,2N+1-j] = 0\) (anti-diagonal entry is zero, since
  \(2N+1-j \neq j\) for all \(j\) in the cell range \(\{1,\ldots,2N\}\)),

so the diagonal solution **violates** (G) by the maximum possible
amount. The relaxation is forced off this point.

This is a *perfectly clean* fix to a problem the project has been
stuck on for four iterations. It works because (G) is a genuine
algebraic identity for the rank-one truth \(Q^* = pp^*\top\) coming
from any *even* admissible \(f\), but the diagonal mixture (which
corresponds to the convex combination \(\mu = \frac1n\sum_l\delta_{e_l}\)
of point-mass \(f\)'s, none of which is even) cannot satisfy it.

### Numerical convergence. ✓

The table

| \(N\) | \(\Omega^{(5)}\) |
|---|---|
| 4 | 1.800 |
| 6 | 1.857 |
| 8 | 1.889 |
| 10 | 1.909 |
| 12 | 1.923 |

shows monotone convergence to 2 from below, matching the
*known* analytic value \(C^{\text{even}}_{6.2} = 2\) achieved by
\(f = 2\cdot\mathbf 1_{[-1/4,1/4]}\) (whose autocorrelation is the
triangle of peak height 2 at \(t=0\)). The empirical convergence is
strong evidence that the LP is a valid relaxation for the
even-restricted constant.

## 4. Issues

### Issue 1 (FATAL — inherited). v5 still inherits the v0 band-constraint hand-wave.

The "Validity of the relaxation" section asserts:

> The v1 constraints (PSD lift, mass, marginals, \(Q\ge 0\), **bands**,
> central window) hold by their respective derivations.

But the v0 band derivation is the same chain
\(\frac12\sum_{j+k\in\{m,m+1\}}p_jp_k = \int_{B_m}(f*f)\,dt \le L\Omega\)
whose middle equality is *false* in general (the LHS is actually
\(\int_{B_m}(\bar f *\bar f)\) for the cell-averaged \(\bar f\)),
flagged as Issue 1 in the v1 verdict and persisting unaddressed
through v2, v3, v4, and now v5.

v5's empirical convergence to the known \(C^{\text{even}}_{6.2} = 2\)
*suggests* that the band constraint is secretly valid (at least for
the even-restricted case), but the verifier's job is to certify the
proof as written, and the band-constraint chain is not certifiable.

**This is the only thing standing between v5 and a clean VALID
verdict.** The (E), (F), (G) additions are perfect. The numerical
results are perfect. The issue is purely inherited.

### Issue 2 (cosmetic, but worth flagging). v5 bounds \(C^{\text{even}}_{6.2}\), not \(C_{6.2}\).

The proposer is *explicitly* honest about this:

> Note that \(C^{\text{even}}_{6.2}\ge C_{6.2}\) (restricting the f's
> gives a larger inf), so a lower bound on \(C^{\text{even}}_{6.2}\)
> does **not automatically** translate into a lower bound on
> \(C_{6.2}\). Nonetheless, v5 demonstrates the relaxation technology
> works in principle and reproduces the known analytic bound
> \(C^{\text{even}}_{6.2}\ge 2\) at the level of the convex program.

The honesty is appreciated and matches the verifier's standard
("don't claim things that may or may not be a valid lower or upper
bound"). v5's status as a *demonstration of technology* rather than a
*direct improvement on the SOTA* is correct.

The current SOTA on \(C_{6.2}\) is 1.28 (per `Problem.md`), and the
analytic value \(C^{\text{even}}_{6.2} = 2\) is much larger, so the
even-restricted bound does not directly improve the project's
numerical target. v6 must extend to general \(f\) for v5's machinery
to matter.

### Issue 3 (algebra error in the v6 plan). The parallelogram identity is off by a factor of 2.

The closing v6-plan section says:
\[
\int f(x)f(-x)\,dx \;=\; \frac12\Big(\int f^2 - \int(f - \check f)^2/2\Big)
\;=\; \frac12\int f^2 - \frac14\int(f-\check f)^2.
\]
But the *correct* parallelogram identity is
\[
\int f(x)f(-x)\,dx \;=\; \int f^2 - \frac12\int(f - \check f)^2.
\]
**Derivation.** \((f-\check f)^2(x) = f(x)^2 - 2f(x)f(-x) + f(-x)^2\),
so \(\int(f-\check f)^2 = 2\int f^2 - 2\int f(x)f(-x)\,dx\), giving
\(\int f(x)f(-x)\,dx = \int f^2 - \tfrac12\int(f-\check f)^2\). The
proof's formula is off by a factor of 2 (it has \(\frac12\int f^2\)
instead of \(\int f^2\)).

This is forward-pointer feedback for v6, not a v5 issue. But the v6
construction will need to use the correct identity.

## 5. Comments

### Comment 1. (G) is the right kind of constraint — it sees \(\Omega\) implicitly via rank-one structure.

In the v3 verdict's discussion of the basis-vector mixture, I noted:

> The fix has to come from a *different kind of constraint* — one that
> mixes the moments with \(\Omega\) in a way the basis-vector mixture
> cannot satisfy.

(G) is exactly this. It's a *quadratic identity* (\(p_j\cdot p_{2N+1-j} = p_j^2\))
that holds for the rank-one truth coming from an *even* admissible
\(f\), but doesn't hold for the convex combination of point masses
that is the basis-vector mixture. The mixture has nonzero diagonal
\(p_j^2 = (1/(2N))^2\) but zero off-diagonal \(p_j p_{2N+1-j} = 0\)
(because the masses sit at different cells). So (G) sees and excludes
the mixture.

The lesson generalizes: any rank-one identity of the form
\(Q[\sigma(j),\sigma(k)] = Q[j,k]\) for some symmetry \(\sigma\) of the
problem will exclude convex combinations that don't respect \(\sigma\).
This is a *symmetry-based* fix to the relaxation. For the even case,
the symmetry is reflection; for the general case, v6 will need a
"symmetry that all admissible \(f\)'s share" — which is harder.

### Comment 2. The v6 plan is on the right track but needs more thought.

The v6 plan is to "adapt the anti-diagonal-equals-diagonal trick to a
constraint that is valid for all admissible \(f\)". The proposed
mechanism:

> by working with both \(f\) and \(\check f(x) = f(-x)\) simultaneously
> and using the fact that
> \(\int f(x)f(-x)\,dx = \frac12(\int f^2 - \int(f - \check f)^2/2)\),
> which couples the evenness defect to a quantity that the relaxation
> can constrain.

is the right idea: you decompose
\(f = f_{\text{even}} + f_{\text{odd}}\) (where
\(f_{\text{even}}(x) = (f(x)+f(-x))/2\) and
\(f_{\text{odd}}(x) = (f(x) - f(-x))/2\)), then bound the
\(f_{\text{odd}}\) defect somehow. **But:** \(f_{\text{odd}}\) is *not*
nonneg in general (it's signed), and \(f = f_{\text{even}} + f_{\text{odd}}\)
combined with \(f \ge 0\) gives \(f_{\text{even}}\ge|f_{\text{odd}}|\).
The constraint \(f_{\text{even}}\ge|f_{\text{odd}}|\) is non-convex
(it's an infinity-norm constraint pointwise), and may be hard to
encode in the LP.

A cleaner direction (suggestion for v6): instead of decomposing \(f\),
directly add the constraint
\(Q[j, 2N+1-j] \le \min(Q[j,j], Q[2N+1-j, 2N+1-j])\) (which is the
weakest form of (G) that holds for *any* nonneg \(f\), via
\(p_j p_{2N+1-j} \le \max(p_j^2, p_{2N+1-j}^2)\) by AM-GM, plus
\(Q\succeq 0\) implications). This is *not* an equality but a
correctly-directed inequality. It rules out the diagonal mixture
(where \(Q[j,2N+1-j] = 0 < \min(Q[j,j], Q[2N+1-j, 2N+1-j]) = 1/(2N)\))
without requiring evenness.

Actually wait — let me re-check that AM-GM step. \(p_j p_{2N+1-j}\le ?\).
AM-GM: \(p_j p_{2N+1-j} \le (p_j^2 + p_{2N+1-j}^2)/2\), which is the
upper bound. The lower bound is just \(p_j p_{2N+1-j}\ge 0\). Neither
of these gives the *equality* constraint of (G).

OK so the AM-GM trick doesn't quite work. v6 will need something more
sophisticated.

**Better suggestion for v6:** introduce the *reflection-pairing
constraint*
\((p_j - p_{2N+1-j})^2 \le \delta\)
for some LP-derived \(\delta\), where \(\delta\) measures the "evenness
defect" of \(f\). For even \(f\), \(\delta = 0\); for general \(f\),
\(\delta\) is bounded by \(\int(f - \check f)^2 = 2\int f^2 - 2\int f\check f\),
which is in turn bounded above by \(2\int f^2\), and \(\int f^2 \le \|f\|_\infty\cdot\|f\|_1 = \|f\|_\infty\),
which is *not* directly controlled by the LP either.

So the v6 path requires more thinking. The honest move might be to
*also* compute \(C^{\text{ren}}_{6.2}\) where "ren" is some other
restriction (e.g., monotone or unimodal) that interpolates between
"even" and "general", and to keep extending. The framework is good;
the algebra needs work.

### Comment 3. Empirical convergence to \(C^{\text{even}}_{6.2} = 2\) is strong evidence.

The five data points in the numerical table show monotone convergence
\(\Omega^{(5)}(N)\to 2\) from below as \(N\) grows. If the LP were
*invalid* (constraints too restrictive, excluding admissible
\(f\)'s), \(\Omega^{(5)}\) could overshoot 2 — but it doesn't. This
is empirical evidence for both:

* The (E), (F), (G) additions correctly characterize the even case
  (no over-restriction).
* The inherited band constraints are *probably* valid in the
  even-restricted case (at least, the LP optimum is consistent with
  it).

The empirical convergence is *not* a proof, and the verdict cannot
treat it as one. But it strongly suggests that an honest proof of the
band constraints is *out there* and just needs to be written. v6
should make a serious attempt to nail down the band proof rigorously
*before* extending to the general case.

## 6. Feedback to the proposer

This is the most substantial iteration in the project so far, and the
(G) construction is genuinely clever. The diagnosis chain
v3 (basis-vector mixture defeats Lasserre) →
v4 (diagonal-\(Q\) degeneracy persists in primal-side relaxation) →
v5 (anti-diagonal-equals-diagonal forces non-degeneracy via even
symmetry) is exactly the right sequence of reasoning.

Three concrete suggestions for v6:

1. **Fix the band constraint *first*.** Five iterations of inheritance
   is enough. The right approach is the one I sketched in the v4
   verdict's Comment 1: introduce cell-supremum variables \(M_j\) with
   \(M_j\ge p_j/L\), then bound \(\int_{B_m}(f*f)\) by a sum over
   \(M_jM_k\)-weighted cell-rectangle areas (an honest *upper* bound
   in the right direction). This makes the band constraint rigorous
   and removes the inheritance issue from every future verdict.
2. **Fix the v6 parallelogram identity.** The correct form is
   \(\int f\check f = \int f^2 - \tfrac12\int(f-\check f)^2\), not
   what the v5 closing paragraph writes.
3. **Don't try to extend (G) directly to general \(f\).** Instead,
   keep (G) as the *even-case* fix, and look for a different symmetry
   that *all* admissible \(f\)'s share. One option:
   \(f \to f\ast \mathbf 1_{[-\epsilon,\epsilon]} / (2\epsilon)\) (a
   small mollification): every admissible \(f\) is the limit of a
   sequence of "approximately translation-invariant" functions, and
   you can encode "approximate translation invariance" as an LP
   constraint at finite \(\epsilon\). The right way to think about
   this is via the *autocorrelation lift*: the autocorrelation
   \(F(t) = \int f(x)f(t-x)\,dx\) is *automatically* even
   (\(F(t) = F(-t)\) since you can substitute \(x\to t-x\)), so
   *every* \(f\) — even or not — yields an even \(F\). Working with
   \(F\) instead of \(f\), the (G) trick becomes a *universal* fact
   about the autocorrelation. **This is the cleanest path forward.**

Suggested v6 structure:

* Lift to an \(F\)-side LP: variables are cell averages \(W_j\) of
  \(F\) on the centered grid \([-1/2, 1/2]\), with \(W_j = W_{2N+1-j}\)
  (since \(F\) is even) and \(\sum W_j = 1\) (since
  \(\int F = (\int f)^2 = 1\)).
* The objective is \(\Omega = \max_j W_j\) (since
  \(\Omega = \|F\|_\infty\)), which is *linear* in \(W\).
* Add the rank-one identification \(W_j = (\hat F\cdot\hat F)\)-style
  constraints from the Fourier side.
* The symmetry \(W_j = W_{2N+1-j}\) is the universal version of (G).

This is the "v4 dual side" plan, executed cleanly. It should give a
real lower bound on \(C_{6.2}\) (not just \(C^{\text{even}}_{6.2}\))
and break the \(\Omega = 1\) trap.

---

**Final verdict: INVALID** (by inheritance only — v5's *new* construction is the most substantial advance in the project so far).
