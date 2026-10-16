# Verdict — Autocorrelation Inequality 6.2, v3

## 1. Summary of the proof

v3 carries forward v0/v1/v2 and adds the standard Lasserre level-2 lift
on the moment variables \(p_j = \int_{I_j} f\):

* **(G) Moment-matrix PSD lift at level 2.** Variables
  \(y_\alpha\) for multi-indices \(|\alpha|\le 4\), with the moment
  matrix \(M_2[i,j] = y_{\alpha^{(i)}+\alpha^{(j)}}\) over multi-indices
  of degree \(\le 2\). Drop rank-one and impose \(M_2\succeq 0\).
* **(H) Localizing matrices for \(p_j \ge 0\).** For each
  \(j \in \{1,\ldots,2N\}\), \(L^{(j)}[i,k] = y_{\alpha^{(i)}+\alpha^{(k)}+e_j}\)
  is PSD because in rank-one truth it equals \(p_j(ww^\top)\) with
  \(p_j\ge 0\).
* **(I) Linear marginals.** \(\sum_j y_{\alpha+e_j} = y_\alpha\) for
  every \(|\alpha|\le 3\), the propagation of \(\sum_j p_j = 1\) to
  higher-order moments.

The band-integral constraints from v0 are inherited.

The proof closes with a **"Numerical observation"** section that
proves a remarkable fact: the LP/SDP at *every* finite Lasserre level
admits a feasible point at \(\Omega = 1\), namely the moment sequence
of the *uniform measure on the basis vectors of the simplex*. So no
Lasserre tightening can move the bound off the trivial value.

## 2. Verdict

**INVALID.**

The v3 additions (G), (H), (I) are *correct* — they are the standard
Lasserre level-2 moment-SOS lift, and the rank-one truth verifies all
of them by direct substitution. **The v3 additions are not the
problem.**

The problem is inherited: v3 carries forward the v0 band-integral
constraint without re-justifying it, and that constraint has been the
unresolved Issue 1 of the verdicts on v1, v2, and now v3. The
proof's chain
\[
\tfrac12\!\!\sum_{j+k\in\{m,m+1\}}\!\! p_jp_k \;=\; \int_{B_m}(f*f)\,dt
\;\le\; L\|f*f\|_\infty
\]
asserts an *equality* in the middle that holds for piecewise-constant
\(f\) only — for general \(f\), the LHS is
\(\int_{B_m}(\bar f*\bar f)\) where \(\bar f\) is the cell-averaged
function, and the two are not equal. The original v1 verdict gave the
explicit counterexample \(f = 4\cdot\mathbf 1_{[-1/8,1/8]}\), \(N=1\),
band \(B_1\): LHS \(= 1/8\), \(\int_{B_1}(f*f) = 0\). Three iterations
later, this gap is still not addressed.

**Important caveat I owe the proposer.** I have tried hard to find an
admissible \(f\) that *actually* violates the LP constraint as a whole
(i.e.\ \(\frac12\sum p_jp_k > L\|f*f\|_\infty\) for some band) and
have failed. So the band constraints may secretly be valid by some
*other* argument the proof doesn't make. But the verifier's job is to
certify the proof as written, not to guess that an unproven step
"should" be true. As written, the band constraint chain is not
rigorous, and so the LP optimum is not certified to be a lower bound
on \(C_{6.2}\). v3 inherits this and is therefore INVALID.

(My counterexample-search difficulty does suggest a path forward for
v4: the band constraint may be salvageable by replacing the false
equality with a *different* upper-bound chain, e.g.\
\(\int_{B_m}(\bar f*\bar f) \le L\|f*f\|_\infty\) proved directly via
a Young-type argument. See §4 below.)

## 3. Verification of the v3 additions (all correct)

### (G) Moment-matrix PSD lift. ✓
The level-2 moment matrix \(M_2[i,j] = y_{\alpha^{(i)}+\alpha^{(j)}}\)
is, in the rank-one truth, the Gram matrix
\((vv^\top)[i,j] = v_iv_j = p^{\alpha^{(i)}+\alpha^{(j)}}\), so PSD.
The relaxation (drop rank-one, keep \(M_2\succeq 0\)) is the standard
Lasserre lift.

### (H) Localizing matrices. ✓
\(L^{(j)} = p_j(ww^\top)\) in the rank-one truth, with \(p_j\ge 0\).
The localizing matrix is the standard SOS dual to the inequality
\(p_j\ge 0\), and produces a real LP constraint
("\(p_j\ge 0\) on the support of the moment measure").

### (I) Linear marginals. ✓
\(\sum_j y_{\alpha+e_j} = \sum_j p_jp^\alpha = p^\alpha\sum_j p_j
= p^\alpha = y_\alpha\), using the mass constraint.

All three constraints hold for the rank-one lift of any admissible
\(f\), so they are valid additions to the relaxation. They cannot
*fail* to be valid — that's the part the proof gets right.

### Inherited issue: the v0 band constraint is still hand-waved.

This is the same issue as v1 Issue 1 and v2 Issue 4. v3 acknowledges
nothing about this and the v0 docstring (which the proof cites)
contains the same false equality. The chain still does not survive
contact with a non-step-function \(f\).

## 4. Comments (especially on the "Numerical observation")

### Comment 1. The basis-vector-measure construction is *correct and important*.

The proof constructs an explicit feasible point of v3's LP at
\(\Omega = 1\):
\[
y_\alpha \;=\; \mathbb E_\mu[p^\alpha], \qquad
\mu = \frac1n\sum_{l=1}^{n}\delta_{e_l},
\]
where \(\mu\) is the uniform measure on the standard basis vectors of
the simplex \(\Delta_n\). I verified this:

* \(y_{\mathbf 0} = 1\), \(y_{e_j} = \mathbb E_\mu[p_j] = 1/n\), so
  \(\sum_j y_{e_j} = 1\) (mass).
* \(y_{e_j+e_k} = \mathbb E_\mu[p_jp_k] = (1/n)\delta_{jk}\) (since
  basis vectors have at most one nonzero coordinate). So
  \(Q[j,k] = (1/n)\delta_{jk}\) is diagonal with trace \(1\).
* The moment matrix \(M_2\) is the rank-\(n\) Gram matrix
  \(\sum_l (1/n) v_l v_l^\top\) where \(v_l\) is the level-\(\le 2\)
  monomial evaluation at \(e_l\). PSD by construction.
* Each localizing matrix \(L^{(j)} = (1/n) v_j v_j^\top\) (where
  \(v_j\) is the level-\(\le 1\) monomial evaluation at \(e_j\)). PSD.
* The marginals (I) hold by construction (they are linear identities
  satisfied by any moment sequence of a probability measure on
  \(\Delta_n\)).
* The band constraint at band \(m\) reduces to
  \((1/2)\sum_{j+k\in\{m,m+1\}}(1/n)\delta_{jk}\). When \(m\) is even
  this is \((1/2)\cdot(1/n) = 1/(2n) = 1/(2 \cdot 2N) = L\); when
  \(m\) is odd, \(0\). In both cases \(\le L\Omega\) with \(\Omega = 1\).

So the basis-vector measure is a *bona fide* feasible point of v3 at
\(\Omega = 1\), and v3's LP optimum is therefore \(\le 1\), the
trivial value. **Confirmed.**

### Comment 2. "No finite Lasserre level can rule this out" is correct.

The proof's structural conclusion is:

> The basis-vector measure \(\mu\) is a *real* probability measure on
> \(\Delta_n\). Therefore its moment sequence \(y_\alpha = \mathbb
> E_\mu[p^\alpha]\) is a valid moment sequence at every level
> \(\ell\) — the level-\(\ell\) moment matrix is the rank-\(n\) Gram
> matrix of the basis-vector evaluations, which is always PSD. No
> finite Lasserre level can rule this out.

This is correct: the Lasserre / moment-SOS hierarchy gradually tightens
the relaxation as \(\ell\) increases, but the *only* way it can exclude
a moment sequence \(y\) is by finding a higher-order PSD constraint
violated by \(y\). Since \(\mu\) is a real measure on the feasible
domain \(\Delta_n\), every moment matrix derived from \(\mu\) at every
level is PSD by definition. So Lasserre tightening *cannot help* here.

This is a project-level insight: it tells you that adding more moment
constraints to the same \((p,Q)\)-side LP cannot break \(\Omega = 1\).
You need a *different* kind of constraint — one that sees \(\Omega\)
itself, not just the moments of \(p\).

### Comment 3. The diagnosis ("must decouple \(\Omega\) from a single \(f\)") is right.

The proof's diagnosis:

> The key fix must therefore decouple \(\Omega\) from a single \(f\),
> e.g. by introducing constraints that involve \(\Omega\) together with
> the moments in a way that no convex combination of point-mass
> distributions can satisfy.

is the right framing. The basis-vector measure is the convex
combination \(\frac1n\sum_l \delta_{e_l}\), and each \(\delta_{e_l}\)
corresponds to "\(f\) is a delta at the centre of cell \(l\)" — for
each such \(f\), the actual \(\|f*f\|_\infty\) is *infinite*, and yet
the convex combination "averages out" to \(\Omega = 1\). The LP's
band constraint is too weak to see this averaging, because it only
sees the second-order moments \(p_jp_k\) and not higher-order
information.

What does work: the **autocorrelation/Fourier side**, which is what
v4 plans. White's Section 4 LP works because the Fourier-side
constraint \(A_{2m}\le 0\) (cosine coefficients of \(M=f*g\)) is a
*linear* condition on the cell averages of \(M\), and is *not*
satisfied by the basis-vector measure (because the basis vectors
correspond to different \(M\)s, and "averaging" them changes
\(A_{2m}\) in a way the LP can detect).

This is the correct direction. v4 should work.

### Comment 4. The inherited band constraint is still the elephant in the room.

Three iterations now have inherited the same hand-waved
band-constraint chain without addressing it. Even if the band
constraint is *secretly* valid, the rigorousproof is supposed to
*prove* it. Some options:

* **Direct \(\bar f \to f\) comparison.** Try to prove
  \(\int_{B_m}(\bar f*\bar f) \le L\|f*f\|_\infty\) directly for any
  admissible \(f\), without going through \(\|\bar f*\bar f\|_\infty\).
  The most promising route is via the Plancherel/Young's chain
  \(\int_{B_m}(\bar f*\bar f) \le \|1_{B_m}\|_2 \cdot \|\bar f*\bar f\|_2
  = \sqrt L \cdot \|\hat{\bar f}\|_4^2\), then bound
  \(\|\hat{\bar f}\|_4^2 \le \|\hat f\|_4^2 \le \sqrt L \cdot \|f*f\|_\infty\)
  (the last step needs work but might go through). I have not verified
  this route.
* **Replace the band integral by a *correctly* discretised one.** Use
  the rectangular cell-area decomposition
  \(K_{m,j,k} = \mathrm{area}\big(I_j\times I_k\,\cap\,\{(x,y):x+y\in B_m\}\big)\),
  and bound
  \(\int_{B_m}(f*f) \le \sum_{j,k}\|f|_{I_j}\|_\infty\,\|f|_{I_k}\|_\infty\,K_{m,j,k}\)
  (which is the *actual* rectangle-by-rectangle upper bound, not the
  step-function equality v0 wrote down). This needs the
  cell-supremum \(\|f|_{I_j}\|_\infty\) as an LP variable, not just
  the cell-average \(p_j\) — i.e.\ a finer discretisation.
* **Drop band constraints entirely** in favour of Fourier-side
  constraints (the v4 plan).

Whatever route v4 takes, **the verdict on the band constraint must be
written down somewhere**. As of v3, it is still in the same hand-waved
state as v0.

## 5. Feedback to the proposer

The Lasserre level-2 additions in v3 are *technically correct and
mathematically clean*. They are a no-op for the bound (because of the
basis-vector measure construction), but as a *framework* addition
they're fine.

The two most important things in this iteration are not the
constraints but the **structural observations**:

1. **Lasserre cannot help here**, because the loose feasible point is
   a real moment sequence and survives every level. This is a
   project-level diagnostic that should make you confident the next
   move is *not* to throw more SDP machinery at the same \((p, Q,
   M_2, \ldots)\) hierarchy.
2. **The fix has to bring \(\Omega\) into the constraints in a
   different way**, e.g., via Fourier-side constraints on \(\hat F(k)\)
   that the basis-vector measure does not satisfy.

Both are right.

For v4, the cleanest move is probably:

1. **Switch to the Fourier-side LP analogous to White's Section 4.**
   Discretize \(F = f*f\) (not \(f\)) by cell averages \(W_j\) of
   \(F\) on \([-1/2, 1/2]\), and use the constraints
   \(\int F = 1\), \(W_j \le \Omega\), and the Fourier-coefficient
   inequalities derived from \(\hat F = \hat f^2\). The
   \(\hat F = \hat f^2\) constraint *is* a real bilinear constraint
   that the basis-vector measure does not satisfy (different \(f\)'s
   give different \(\hat f^2\)).
2. **Resolve the band constraint properly**, ideally by switching to
   the cell-average \(W_j\) of \(F\) directly (so the band integrals
   are *exact*, not approximated through cell averages of \(f\)).
3. **Avoid the v2 (F) Gibbs trap.** Use a Fejér kernel or the
   coefficient-side formulation; do not test partial Fourier sums
   pointwise.

If v4 does this, the band-constraint hand-wave goes away (replaced by
the exact \(W_j\) cell averages of \(F\)) and the basis-vector
degeneracy is broken (because \(\hat F = \hat f^2\) is rank-one in
\((a_k, b_k)\) and the basis-vector measure does not lift to a rank-one
sequence in those variables).

---

**Final verdict: INVALID.**
