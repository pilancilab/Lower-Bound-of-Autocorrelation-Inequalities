# Verdict — Autocorrelation Inequality 6.2, v4

## 1. Summary of the proof

v4 inherits the v0/v1 \((p, Q)\) SDP lift and adds a new family of
constraints (W) — *windowed integral bounds* — derived as follows:

For any window \(B = [a, b]\subset[-1/2,1/2]\) with full width \(h_B = b - a\),
\[
\int_B (f*f)(t)\,dt = \iint f(\tau)f(x)\,\mathbf 1_{a\le\tau+x\le b}\,d\tau\,dx
\;\le\; h_B\cdot\Omega,
\]
and the iterated integral is bounded *from below* by the sum over
fully-contained cells:
\[
\iint f(\tau)f(x)\,\mathbf 1_{a\le\tau+x\le b}\,d\tau\,dx
\;\ge\; \sum_{(j,k):\,I_j\times I_k\subset\text{preimage}} p_j p_k.
\]
Combining gives the LP constraint
\[
\sum_{(j,k):\,I_j\times I_k\subset\text{preimage}} Q[j,k] \;\le\; h_B\,\Omega.
\tag{W}
\]
Windows used: centered \([-h, h]\) and shifted \([\pm 1/4 - h, \pm 1/4 + h]\)
for half-widths \(h\in\{2L, 3L, \ldots, 1/2\}\).

The proof closes with a "Numerical observation": even with (W), MOSEK
returns \(\Omega^{(4)} = 1\), the trivial bound; the diagonal \(Q\)
solution evades the (W) constraints exactly as it evaded v1/v2/v3.

The proof's stated v5 plan is to *exhibit a dual feasible point* — a
positive measure \(\nu\) on \([-1/2, 1/2]\) — for which
\(\inf_f\int(f*f)\,d\nu/\int d\nu > 1\). This is a structural change of
direction.

## 2. Verdict

**INVALID — by inheritance, not by anything new.**

The v4 (W) constraint is **mathematically valid** and is, in fact, the
**first time in the project** that a band-style constraint has been
proved without the false-equality hand-wave that infected v0/v1/v2/v3
(see Comment 1 below). The v4 *additions* are rigorous.

However, the proof's "Validity of the relaxation" section explicitly
re-asserts that v4 inherits all v0 constraints, including the hand-waved
band constraints from v0. So v4's claim
\(\Omega^{(4)}\le C_{6.2}\) still depends on the unresolved v0 band
constraint, and is therefore not certifiable as written.

This is the same inheritance issue as v1, v2, v3. v4 is **closer to a
fix** than any prior iteration (because (W) is a *correct* alternative
to the v0 bands), but does not actually drop the inherited bands.

## 3. Verification

### Constraint (W) is rigorously derived. ✓

The chain
\[
\int_B(f*f)\,dt \;=\; \iint f(\tau)f(x)\,\mathbf 1_{a\le\tau+x\le b}\,d\tau\,dx \;\ge\; \sum_{\text{full}} p_j p_k
\]
is correct, because:
* The first equality is Fubini on the integrand
  \(\int f(\tau) f(t-\tau)\,d\tau\) (substituting \(x = t - \tau\)).
* The second inequality follows because for a *fully-contained* cell,
  the indicator is identically 1, so
  \(\int_{I_j\times I_k} f(\tau)f(x)\cdot 1\,d\tau\,dx = (\int_{I_j}f)(\int_{I_k}f) = p_j p_k\)
  *exactly*. Partially-contained and excluded cells contribute
  nonneg amounts to the iterated integral but are not in the sum, so
  dropping them gives a lower bound.
* The upper bound \(\int_B(f*f)\le h_B\cdot\|f*f\|_\infty\) is the
  trivial \(\int 1_B\cdot F \le \|1_B\|_1\|F\|_\infty\).

So the (W) constraint at any window \(B\) holds for every admissible
\(f\) and the rank-one lift \(Q^* = pp^*\top\). ✓

### What's *new* about (W) compared to v0's band constraints.

The crucial difference: v0's band constraint at a band of width \(L\)
asserts
\[
\tfrac12\!\!\sum_{j+k\in\{m,m+1\}}\!\! Q[j,k] \;\le\; L\,\Omega
\]
which the v0 file derives via a chain that includes the false equality
\(\frac12\sum p_jp_k = \int_{B_m}(f*f)\). The factor \(1/2\) comes from
*partial* triangular intersections of the cell rectangles \(I_j\times I_k\)
with the band, which the v0 derivation handled by an unjustified
"step-function" pretence.

v4's (W) at the same width \(L\) avoids this entirely by *only counting
fully-contained cells* — but at width \(L\), no cell is fully contained
(\(\inf(\tau+x) = -1/2 + (j+k-2)L\), \(\sup(\tau+x) = -1/2 + (j+k)L\),
spanning an interval of length \(2L\) which cannot fit in a width-\(L\)
window). So **(W) at width \(L\) is the trivial \(0\le L\Omega\)** —
rigorous but useless. (W) only becomes non-trivial at widths \(\ge 2L\),
where one or more \(j+k\) values become fully contained.

So (W) and v0's bands are *not comparable*: (W) is rigorous but
trivially-weak at fine scales; v0's bands are stronger at fine scales
but unproven.

**The right move (which v4 doesn't make):** drop the v0 bands entirely
and replace them by (W) at multiple widths. The resulting LP would be
*rigorously* a lower bound on \(C_{6.2}\), at the cost of being weaker
(possibly losing the trivial \(\Omega\ge 1\) — needs checking). v4
keeps both, which is INVALID by inheritance even though (W) on its
own would be VALID.

### Numerical observation: small arithmetic error, correct conclusion.

I checked the diagonal-\(Q\) feasibility argument carefully. With
\(p = (1/(2N))\mathbf 1\), \(Q = (1/(2N))I\), and a centered window
\([-h, h]\) of full width \(2h\):

* Diagonal cells fully contained: those \(j\) with
  \(I_j\subset [-h/2, h/2]\). Number is approximately
  \(h/L = 4Nh\), not \(2hN\) as the proof says. (The proof appears to
  use \(L = 1/(2N)\), which would correspond to \(N\) cells, not
  \(2N\); the v0 setup has \(L = 1/(4N)\) and \(2N\) cells.)
* Sum of diagonal entries on these cells:
  \((4Nh)\cdot(1/(2N)) = 2h\).
* Constraint becomes \(2h \le 2h\cdot\Omega\), i.e., \(\Omega\ge 1\).

The proof says the constraint becomes "\(h_B \le 2 h_B \Omega\), i.e.,
\(\Omega\ge 1/2\)". With \(h_B\) reinterpreted as half-width, the proof
has \(\sum = h_B\) but I get \(\sum = 2h_B\). **Off by a factor of 2.**

This is a cosmetic error: the proof's *conclusion* ("(W) doesn't beat
\(\Omega = 1\) on the diagonal solution") is correct — both my \(\Omega\ge 1\)
and the proof's \(\Omega\ge 1/2\) are below the binding band-constraint
value of 1. But the specific constant should be \(\Omega\ge 1\), not
\(\Omega\ge 1/2\).

### The structural conclusion is correct.

The proof's deep diagnosis:

> linear and PSD relaxations on \((p, Q)\) alone cannot exclude the
> basis-vector mixture, and the relaxation cannot strictly improve on
> \(\Omega = 1\) without introducing constraints that fundamentally
> distinguish a single \(f\) from a mixture of delta-supported \(f\)'s.

is **correct and is the same insight as v3's "Numerical observation"
about Lasserre levels.** The basis-vector measure
\(\mu = \frac1n\sum_l\delta_{e_l}\) on \(\Delta_n\) gives a *real*
moment sequence at every level, and any constraint expressible purely
in terms of moments \(y_\alpha\) is automatically satisfied. (W) only
involves the second moments \(p_jp_k\), so it falls into this category.

The fix has to come from a *different kind of constraint* — one that
mixes the moments with \(\Omega\) in a way the basis-vector mixture
cannot satisfy.

## 4. Comments

### Comment 1. v4's (W) is *almost* the right way to fix v0's bands.

v4 has discovered the *correct* derivation pattern:
\(\int_B(f*f) = \iint f f \mathbf 1 \ge \sum_{\text{full}} p_jp_k\),
followed by \(\int_B(f*f) \le h_B\Omega\). This is the rigorous,
hand-wave-free version of the band constraint, and it has been missing
since v0.

The main thing v4 lacks is to *also* bound the partial-cell
contributions from above, so as to recover a v0-strength constraint at
fine scales. The standard way to do this is:
\[
\int_B(f*f) \;\le\; \sum_{(j,k):\,I_j\times I_k\,\cap\,\text{preimage}\neq\emptyset} \big(\sup_{(τ,x)\in I_j\times I_k}f(\tau)f(x)\big)\cdot \mathrm{area}\big(I_j\times I_k\,\cap\,\text{preimage}\big),
\]
which gives a *correct upper bound* in terms of cell suprema (not just
cell averages). The cell suprema are themselves not LP variables, but
they can be bounded by introducing additional variables
\(M_j := \sup_{x\in I_j} f(x)\) with linear constraints
\(p_j = \int_{I_j}f \le M_j L\) (i.e., \(M_j \ge p_j/L\), the
cell-average lower bound). Plugging back in gives a constraint on
\(p_j p_k\) and \(M_j M_k\) that is rigorous *and* uses the
partial-cell areas.

This is more work than v4 puts in, but it's the missing ingredient
that would let v4 *replace* the v0 bands rigorously instead of
inheriting them.

### Comment 2. The v5 dual-side plan is the right structural move.

The proof's stated v5 plan:

> v5 will attempt this distinction by working on the dual side: rather
> than relax the primal, we exhibit a dual feasible point (a positive
> measure \(\nu\) on \([-1/2, 1/2]\)) for which
> \(\inf_f \int(f*f)\,d\nu/\int d\nu > 1\), giving an honest lower bound.

is **exactly right**. This parallels what 6.3 v3 needs to do (switch
to a concave-objective formulation) and what White does in his
Section 4 LP for the Erdős problem (the LP is built by *choosing*
linear functionals of \(M\) that admit Fourier control, not by
relaxing the inf-sup directly).

Concretely: for any positive measure \(\nu\) on \([-1/2, 1/2]\),
\[
\int(f*f)\,d\nu \;=\; \iint f(\tau)f(x)\,\nu\big([\,(\tau+x):\,\cdot\,]\big)\,d\tau\,dx
\]
and you want to bound this *from below* by a function of the moments
of \(f\), then divide by \(\int d\nu\). If you can find any \(\nu\)
for which the resulting infimum exceeds 1, you have proven
\(C_{6.2} > 1\).

A natural starting point: take \(\nu(t) = \cos^2(\pi t)\,dt\) (or
similar nonneg trig polynomial) and use Plancherel to express the
inner-product side as \(\sum_k |\hat\nu(k)|^2 |\hat f(k)|^4\), which
is a linear functional in the *Fourier coefficients* of \(f\) — the
exact same variables you'd use in a Fourier-side relaxation. The
basis-vector mixture *does not* satisfy these Fourier constraints,
because the Fourier coefficients of a delta-supported \(f\) are
explicit oscillating functions of the delta location, not constants.

So v5 should break the bound. Good plan.

### Comment 3. The proof's arithmetic in the Numerical observation should be fixed.

The line "\(\sum_{D_B} Q[j,j] = 2 h_B \cdot N \cdot 1/(2N) = h_B\)"
should be "\(= 2h_B\)" (and the conclusion "\(\Omega\ge 1/2\)" should
be "\(\Omega\ge 1\)"). The derivation uses inconsistent
half-width-vs-full-width conventions for \(h_B\). The conclusion
("(W) is consistent with \(\Omega = 1\)") survives the fix; it just
becomes "(W) gives the same \(\Omega\ge 1\) as the v0 bands, *not*
strictly weaker".

## 5. Feedback to the proposer

This is real progress. Two concrete observations:

1. **You have re-derived the band constraint correctly for the first
   time.** The chain
   \(\int_B F = \iint f f\mathbf 1 \ge \sum_{\text{full}} p_jp_k\)
   followed by \(\int_B F\le h_B\Omega\) is the *right* way to do
   this. It's been wrong in v0/v1/v2/v3, and you've fixed it. Now
   *use* the fix: drop the v0 bands and replace them with (W). The
   resulting LP will be rigorous (no inherited hand-wave) at the cost
   of being weaker at fine scales. To recover the strength, add cell
   suprema \(M_j\) with the bound \(M_j \ge p_j/L\) and use them to
   handle partial cells (Comment 1).
2. **The v5 dual-side plan is correct and is also the right move for
   problem 6.3.** Both problems get stuck on the same structural
   issue: linear-in-moments constraints alone cannot break the
   basis-vector / convex-combination loose feasible point. Both
   problems need a *concave* (or dual-side) reformulation of the
   objective to make this work. The dual-side angle you're proposing
   for v5 is the natural counterpart to "switch the objective to a
   concave surrogate" that I've been suggesting for 6.3.

A concrete suggestion for v5:

* Pick \(\nu(t) = \frac{1}{Z}(1 - 4|t|)_+\) (a triangular kernel
  centered at 0, supported on \([-1/2, 1/2]\), with \(Z = 1/2\) so
  \(\int d\nu = 1\)). Then
  \(\int(f*f)\,d\nu = \frac{1}{Z}\iint f(\tau)f(x)(1 - 4|\tau+x|)_+\,d\tau\,dx\),
  which is *bilinear* in \(f\) and admits the same SDP lift you've
  been using.
* Compute the value of this functional on the rank-one lift
  \(Q^* = pp^\top\). It will be a quadratic form in \(p\) with a
  *concave* dependence on \(p\) (because the kernel \((1 - 4|t|)_+\)
  is a positive measure, hence its Fourier transform is positive,
  and the resulting inner product is sum-of-squares-like).
* Now the relaxation maximizes this quadratic functional over
  \((p, Q)\) subject to \(Q\succeq pp^\top\), \(Q\ge 0\), \(\sum p = 1\),
  and the (W) constraints. This is a SDP, and the basis-vector
  mixture does not maximize it — because the diagonal \(Q\) gives a
  smaller value than a non-diagonal one.
* If the SDP optimum is \(\le c\), then \(\int F\,d\nu \le c\) for
  every admissible \(f\), hence \(\Omega = \|F\|_\infty \ge \int F\,d\nu / Z = c \cdot 2\). If you can choose \(\nu\) so that this gives
  \(c \cdot 2 > 1\), you have a lower bound \(C_{6.2} > 1\).

This is a viable v5. Worth trying.

---

**Final verdict: INVALID** (by inheritance only — v4's *new* constraint is the first rigorous band-style constraint in the project).
