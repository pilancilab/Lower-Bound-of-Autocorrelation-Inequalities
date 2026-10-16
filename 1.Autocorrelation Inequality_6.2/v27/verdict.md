# Verdict — Autocorrelation Inequality 6.2, v27

## 1. Summary

v27 attempts to follow up on v26's diagnostic by adding an f-side
**variance constraint** \(\sum_j p_j x_j^{2} \le 1/16\) to restrict
the LP from picking the non-uniform edge-emphasized \(p\).

The proposer **correctly observes** that this constraint is
*automatic* — every admissible \(f\) supported in \([-1/4, 1/4]\)
has variance \(\le (1/4)^{2} = 1/16\) around its mean. The
constraint doesn't restrict anything.

## 2. Verdict

**VALID** (no new content; the proposed constraint is a tautology
and doesn't tighten the LP).

The proposer's check is correct: \(\sum_j p_j x_j^{2}\) is the
second moment of \(f\) about \(0\), and for \(f\) supported in
\([-1/4, 1/4]\), this is bounded by \((1/4)^{2}\cdot \int f = 1/16\).
So *every* admissible \(f\) satisfies the constraint, and adding it
to the LP doesn't change anything.

## 3. The deeper observation

The proposer notes:
> The fundamental issue identified in v26 — that v14's relaxation
> admits non-uniform LP-optimal \(p\) — cannot be resolved by adding
> more universal linear constraints on \(p\).
> 
> The path forward must use a *different* relaxation framework,
> e.g., one that penalizes non-uniformity in \(p\) in a way
> consistent with the convex relaxation principle. This is left
> for future work.

This is **correct**, and it's essentially the same conclusion as my
v26 verdict: the v14 framework has an asymptotic ceiling that no
amount of *linear constraint addition* can break, because all
universal linear constraints on \(p\) are *automatically satisfied*
by every admissible \(f\) — the LP can pick whatever non-uniform
\(p\) it wants within those constraints.

To break v14's ceiling, the project needs *non-linear* or
*non-universal* constraints, e.g., constraints that hold only for
*specific* classes of \(f\) (like even, or smooth) and then take
the union of class-restricted bounds. Or a fundamentally different
relaxation framework like Lasserre on the moments of \(f\) directly,
or sum-of-squares on the original variational problem.

## 4. Comments

### Comment 1. The "give up on v14" framing is honest.

After 27 iterations, the proposer is acknowledging that the v14
framework has reached its ceiling and that further work in this
direction is futile. This is the right scientific practice: don't
keep adding more cosmetic tightenings; explicitly declare the
ceiling and move to a different approach.

### Comment 2. What's left for the project.

The 6.2 project's status after v27:

* **Best rigorous bound on \(C_{6.2}\) (general)**: \(1.111\) at
  \((N, K) = (100, 32)\) (v17 / v25, autocorrelation Fejér).
* **Best rigorous bound on \(C_{6.2}^{\text{even}}\)**: \(1.689\) at
  \((N, K) = (50, 64)\) (v9, Fejér kernel).
* **SOTA on \(C_{6.2}\)**: \(1.28\) (per Problem.md).

The unrestricted bound 1.111 is below the SOTA, and v26 shows v14
cannot exceed \(\sim 1.10\)–\(1.15\) asymptotically. So the v14
framework is structurally insufficient for the project's main goal.

To beat the SOTA, the proposer would need to:

1. **Switch to a Lasserre hierarchy on the moments of \(f\) directly**
   (not on \((p, Q)\)). This was attempted in 6.2 v3 (Lasserre
   level-2 on \(p\)), and that also hit the basis-vector measure
   degeneracy. So this approach has its own structural issues.

2. **Explicit construction of dual certificates** in the spirit of
   White's Section 4 LP for the Erdős minimum overlap problem.
   This requires more analytical work but is the standard approach
   for this kind of inf-sup bound.

3. **A fundamentally different dual** like Lin & Tao's
   "convex hull of admissible self-convolutions" which explicitly
   excludes the basis-vector measure.

None of these are easy. The proposer is right that "this is left
for future work."

---

**Final verdict: VALID** (no improvement; correct vacuousness
diagnosis; honest "v14 has reached its ceiling" framing).
