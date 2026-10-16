# Verdict — Autocorrelation Inequality 6.2, v1

## 1. Summary of the proof

The proof works in the normalization
\[
f\ge 0,\qquad \operatorname{supp}(f)\subset[-1/4,1/4],\qquad \int f = 1,
\]
discretizes \([-1/4,1/4]\) into \(2N\) equal cells \(I_j\) of width \(L=1/(4N)\),
and introduces the lifting variables
\[
p_j := \int_{I_j} f \;\ge\; 0,
\qquad Q[j,k] := p_j p_k.
\]
It then claims that two new linear constraints can be added to the v0 SDP
relaxation while preserving the property that the relaxed optimum
\(\Omega^{(1)}\) is a lower bound on \(C_{6.2}\):

* **(A) Element-wise nonnegativity:** \(Q[j,k]\ge 0\).
* **(B) Marginal/row-sum:** \(\sum_{k=1}^{2N} Q[j,k] = p_j\).

Both are proved by direct substitution of \(Q[j,k]=p_j p_k\) and the
mass condition \(\sum_k p_k = 1\). The proof closes by exhibiting the
rank-one lift \((p^\ast,Q^\ast)\) coming from any admissible \(f\),
asserting that it is feasible for the full v1 program, and concluding
\(\Omega^{(1)}\le \|f*f\|_\infty\), hence \(\Omega^{(1)}\le C_{6.2}\).

## 2. Verdict

**INVALID.**

Constraints (A) and (B), in isolation, are correctly proven for the
rank-one lifting and would be valid tightenings *of an already-valid
program*. The proof, however, also makes a global validity claim
("the v1 optimum is a valid lower bound on \(C_{6.2}\)") whose
justification depends on a step that is mathematically wrong as
written. Details below.

## 3. Issues

### Issue 1 (load-bearing). The "band-integral identity" is false for general \(f\).

The closing argument writes
\[
\frac{1}{2}\!\!\sum_{j+k\in\{m,m+1\}}\!\! Q^\ast[j,k]
\;=\; \int_{B_m}(f*f)(t)\,dt
\;\le\; |B_m|\,\|f*f\|_\infty
\;=\; L\,\Omega^\ast,
\]
with the parenthetical that the equality "relies only on \(f\) being a
step function on the grid; in the continuous case the equality is
replaced by an inequality after standard piecewise approximation."

This is exactly the kind of step the verifier instructions forbid. Let
me make it precise.

Define the cell-averaged function
\[
\bar f(x) := \frac{p_j}{L}\quad\text{for }x\in I_j
\qquad (\bar f \text{ is the conditional expectation of } f
   \text{ given the partition}).
\]
A direct computation (split each rectangle \(I_j\times I_k\) into the
two triangles cut out by \(t\in B_m\)) gives the *correct* identity
\[
\frac{1}{2}\sum_{j+k\in\{m,m+1\}} p_j p_k
\;=\; \int_{B_m}(\bar f*\bar f)(t)\,dt,
\]
**not** \(\int_{B_m}(f*f)(t)\,dt\). So the LHS of the proof's chain is
the band integral of \(\bar f*\bar f\), not of \(f*f\). The proof has
silently replaced \(f\) by its piecewise-constant projection. The two
quantities are not equal in general, and — crucially — they are *not
even comparable in a fixed direction*.

**Concrete counterexample to the "equality replaced by inequality" claim.**
Take \(N=1\), \(L=1/4\), so two cells \(I_1=[-1/4,0]\), \(I_2=[0,1/4]\).
Let
\[
f(x) = 4\cdot \mathbf 1_{[-1/8,\,1/8]}(x).
\]
Then \(f\ge 0\), \(\operatorname{supp}(f)\subset[-1/4,1/4]\), \(\int f = 1\), so
\(f\) is admissible. We have \(p_1=p_2=1/2\), and:

* \(f*f\) is a triangle of base \([-1/4,1/4]\) and peak height \(4\). In
  particular, \(f*f\equiv 0\) on \(B_1=[-1/2,-1/4]\), so
  \(\int_{B_1}(f*f) = 0.\)
* \(\bar f \equiv 2\) on \([-1/4,1/4]\), so \(\bar f * \bar f\) is a
  triangle of base \([-1/2,1/2]\) and peak height \(2\). Direct
  integration:
  \[
  \int_{B_1}(\bar f*\bar f) = \int_{-1/2}^{-1/4}2(1+2t)\,dt = \tfrac{1}{8}.
  \]
* And the proof's LHS is
  \(\tfrac{1}{2}\sum_{j+k\in\{1,2\}} p_j p_k = \tfrac{1}{2}(p_1 p_1) = \tfrac{1}{8}\),
  matching \(\int_{B_1}(\bar f*\bar f)\) and **not** matching
  \(\int_{B_1}(f*f) = 0\).

So the proof's claimed equality fails, and the inequality direction
runs the *wrong* way for this band: the LHS of the proof's chain
exceeds \(\int_{B_m}(f*f)\). The "standard piecewise approximation"
appeal cannot rescue this — it asserts an inequality direction that
this example explicitly violates.

(The v1 program's *outer* inequality
\(\frac{1}{2}\sum p_j p_k \le L\,\|f*f\|_\infty\) may still be true for
this particular example because the right-hand side is large enough
to absorb the gap. That is irrelevant: it would have to be proved by
a different argument, and the proof gives no such argument. As
written, the chain is broken by the false middle step.)

### Issue 2. There is no v0.

The proof and the code both reference a v0 baseline ("In v0 the only
constraints on \(Q\) were …"), but
`Autocorrelation Inequality_6.2/` contains only `Problem.md`,
`proposal.md`, and `v1/`. There is no `v0/` folder, no
`v0/convex_code.py`, and no `v0/rigorousproof.md`. So the band
constraints cannot be "inherited from v0" — they must be re-justified
inside `v1/rigorousproof.md`, and they are not. This compounds
Issue 1: a missing baseline plus a hand-waved bridge from the
discrete sum to the continuum integral leaves the entire validity
claim unsupported.

### Issue 3. Two implicit reductions are not justified.

The proof opens with "We work in the normalized setting \(f\ge 0\),
\(\operatorname{supp}(f)\subset[-1/4,1/4]\), \(\int f=1\)," but neither reduction
is derived from the original problem statement (which only requires
\(f\ge 0\) on \(\mathbb R\)). Both are in fact true:

* Replacing \(f\) by \(f\cdot \mathbf 1_{[-1/4,1/4]}\) leaves
  \(\int_{-1/4}^{1/4} f\) unchanged and pointwise *decreases*
  \(f*f\) (convolution of nonnegatives is monotone), hence cannot
  increase \(\max_{|t|\le 1/2}(f*f)(t)\); so the support reduction is
  WLOG.
* The ratio is \(0\)-homogeneous in \(f\), so we may rescale to
  \(\int f=1\).

These are short, but a "rigorous proof" must include them. As
written, a reader would have to take them on faith.

### Issue 4 (minor). (B) is stated as "linear" but is non-trivially restrictive.

It is worth being explicit that combined with the PSD lift and
\(\sum p_k = 1\), constraint (B) forces
\(\mathbf 1^\top Q\,\mathbf 1 = (\mathbf 1^\top p)^2 = 1\), i.e. the all-ones
vector lies in the kernel of \(Q - pp^\top\). The proof presents (B)
as if it were a routine consequence of \(Q=pp^\top\) (which it is for
the rank-one lift), without noting that for the relaxation it is
genuinely tightening. This is not an error, just a missed point that
weakens the case for adding the constraint.

## 4. What would make this VALID

The proposal can be salvaged with the following changes.

1. **Either prove the band constraints from scratch** for general
   \(f\) (not just step functions), or **explicitly cite a v0 proof**
   that does so — and produce that v0 file. In particular, the
   correct route is to bound \(\int_{B_m}(f*f)\) by a cell-wise sum
   *with the right inequality direction* (e.g., by majorizing the
   triangular kernel within each \(I_j\times I_k\) rectangle by its
   maximum on the rectangle, which gives an *upper* bound on the
   true band integral and is then majorized by \(L\|f*f\|_\infty\)).
   Or use a different right-hand side: \(\int_{B_m}(\bar f*\bar f)\)
   only matches \(\int_{B_m}(f*f)\) when \(f\) is already a step
   function on the grid.

2. **Justify the two reductions** (compactification of support and
   mass normalization) at the top of the proof. Two short paragraphs
   suffice; both reductions are correct.

3. **Restate the conclusion modestly:** "Constraints (A) and (B) hold
   for the rank-one lift coming from any admissible \(f\), so adding
   them to a *previously valid* program preserves validity." Then
   the v1 file's job is to certify (A) and (B), and the validity of
   the underlying band program is the v0 file's job. Right now the
   v1 proof tries to do both and fumbles the second.

4. **Mention the kernel-of-\(Q-pp^\top\) consequence of (B)** so the
   reader can see that (B) is genuinely cutting feasible space and
   not merely rephrasing the PSD lift.

## 5. Feedback to the proposer

The proposer is on the right track: lifting to a moment matrix \(Q\)
and adding necessary linear/PSD conditions on \(Q\) is exactly the
template White uses (Lemmas 3, 5, 7 in 2201.05704v1). The proposer
should note, however, that White is *very* careful about two things
that are missing here:

* **Tail / approximation control is always quantitative.** Lemma 4 of
  White does not say "the tail is small for large \(T\)"; it gives an
  explicit constant \(\frac{1}{4-m^2/T^2}\cdot \frac{2m}{\pi\sqrt{6T^3}}\).
  Any time you replace a continuum object by a finite-dimensional
  one, you owe a bound of the same flavor. The "standard piecewise
  approximation" sentence in this proof is the symptom of having
  skipped that step.
* **The relaxation direction must be checked at every step.** Section
  5 of White lists every constraint and says exactly which lemma
  produces it and why it has the correct direction. A `Constraint
  feasibility` table at the end of `rigorousproof.md` (the analog of
  White's Table 1) would force the proposer to confront exactly this
  question, and would have caught Issue 1 immediately.

Suggested next iteration: split this into v1a and v1b. v1a only adds
constraint (A) (which is bullet-proof). v1b only adds constraint (B)
(also fine, but document the kernel consequence). And, separately,
fix the v0 baseline so that the band constraint is actually proved.

---

**Final verdict: INVALID.**
