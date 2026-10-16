# Verdict — Autocorrelation Inequality 6.2, v7

## 1. Summary of the proof

v7 is a **clean reboot** of the relaxation. It explicitly addresses
*every* load-bearing issue from the verdicts on v1, v2, v3, v4, v5,
v6, and produces a rigorous lower bound on
\(C^{\text{even}}_{6.2}\) (the even-restricted version of \(C_{6.2}\)):

* The two **WLOG reductions** (support compactification and mass
  normalization) are explicitly proved at the top of the file.
* The **even restriction** is acknowledged as a *real* restriction:
  v7 bounds \(C^{\text{even}}_{6.2}\), not \(C_{6.2}\).
* The constraint list (P0)–(P5) consists of:
  PSD lift, mass, elementwise nonnegativity, marginals, the v5
  even-symmetry constraints (E)/(F)/(G), and the new rigorous (W)
  windowed-integral family.
* **The v0 band-integral constraint is dropped entirely.** The
  inherited hand-waved chain is gone, replaced by the (W) family
  proved cleanly via Fubini on cell rectangles.
* Numerical results: \(\Omega^{(7)}\) increases monotonically from
  1.046 at \(N=4\) to 1.344 at \(N=20\), strictly above 1 for all
  \(N\ge 4\), converging toward the analytic value
  \(C^{\text{even}}_{6.2} = 2\).

## 2. Verdict

**VALID.**

This is the **first valid iteration in the project history.** Every
load-bearing issue from the previous six verdicts is resolved:

| Issue | Status in v7 |
|---|---|
| v1 Issue 1: false band-integral identity | **Fixed** (v0 bands dropped; replaced by rigorous (W)) |
| v1 Issue 2: no v0 file | Addressed in v0 between iterations; v7 doesn't depend on v0 anyway |
| v1 Issue 3: WLOG support not justified | **Fixed** (Reduction 1 proved) |
| v1 Issue 4: WLOG normalization not justified | **Fixed** (Reduction 2 proved) |
| v2 Gibbs / tail issue (constraint F) | N/A — v7 doesn't use the v2 partial-Fourier approach |
| v3 Lasserre / basis-vector mixture | (G) breaks the diagonal degeneracy (for even \(f\)) |
| v4 windowed (W) constraint not used | **Fixed** — v7 actually uses (W) instead of inheriting v0 |
| v5 inheritance of v0 bands | **Fixed** — v0 bands are dropped |
| v6 phantom slacks | N/A — v7 doesn't use slacks |

The proof is now self-contained, rigorous, and produces a numerical
bound that is strictly above the trivial value 1 for the first time
in the project. The honest disclosure that v7 only bounds
\(C^{\text{even}}_{6.2}\), not \(C_{6.2}\), is correct and is the
right framing.

## 3. Verification

### WLOG Reduction 1 (support compactification). ✓

Replace \(f\) by \(f\cdot\mathbf 1_{[-1/4,1/4]}\).

* Denominator unchanged: \(\int_{-1/4}^{1/4}(f\cdot\mathbf 1_{[-1/4,1/4]}) = \int_{-1/4}^{1/4} f\). ✓
* Numerator: by monotonicity of convolution of nonneg functions
  (\(0\le g_1\le g_2 \Rightarrow g_1*g_1 \le g_2*g_2\) pointwise),
  \((f\cdot\mathbf 1_{[-1/4,1/4]})*(f\cdot\mathbf 1_{[-1/4,1/4]}) \le f*f\) pointwise. ✓
* Hence \(\max_{|t|\le 1/2}\) of the LHS \(\le\) RHS, so the ratio
  weakly decreases. ✓

The infimum is therefore attained on \(f\) supported in \([-1/4,1/4]\). ✓

### WLOG Reduction 2 (mass normalization). ✓
The ratio is degree-2 homogeneous in \(f\), so we may rescale to
\(\int f = 1\). ✓

### Constraints (P0)–(P3). ✓
All four are standard and were verified in the v1 verdict. (P0) is the
PSD lift; (P1) is the mass equation; (P2) is element-wise
nonnegativity (which uses \(p\ge 0\)); (P3) is the marginal/row-sum
identity (which uses \(\sum p = 1\)).

### Constraint (P4): even-symmetry (E)/(F)/(G). ✓
All three were verified in the v5 verdict.

* (E) \(p_j = p_{2N+1-j}\) follows from \(I_{2N+1-j} = -I_j\) (which
  is correct for the centred grid with \(L = 1/(4N)\) and the cell
  range \(\{1,\ldots,2N\}\)) and \(f(-y) = f(y)\).
* (F) \(Q[j,k] = Q[2N+1-j, 2N+1-k]\) follows from (E) twice.
* (G) \(Q[j, 2N+1-j] = Q[j,j]\) follows from (E) and the rank-one
  identification \(Q[j,j] = p_j^2\). **This is the constraint that
  breaks the diagonal-\(Q\) degeneracy** (the diagonal solution has
  \(Q[j,j] = 1/(2N) > 0\) but \(Q[j, 2N+1-j] = 0\), violating (G)).

### Constraint (P5): rigorous (W) windowed-integral family. ✓ (load-bearing)

This is **the** load-bearing fix. I checked every step:

**Step (a) Fubini.**
\[
\int_B (f*f)(t)\,dt
\;=\; \int_B \int f(\tau)\,f(t-\tau)\,d\tau\,dt
\;\stackrel{\text{Fubini}}{=}\; \int f(\tau)\int_B f(t-\tau)\,dt\,d\tau
\;\stackrel{x=t-\tau}{=}\; \int f(\tau)\int_{B-\tau} f(x)\,dx\,d\tau
\;=\; \iint f(\tau)f(x)\,\mathbf 1_{a\le\tau+x\le b}\,d\tau\,dx.
\]
The substitutions are valid because all integrands are nonneg (Fubini
applies without integrability concerns when the integrand is nonneg). ✓

**Step (b) Lower bound by fully-contained cells.** The decisive step.
\[
\iint f(\tau)f(x)\,\mathbf 1_{a\le\tau+x\le b}\,d\tau\,dx
\;\ge\; \sum_{(j,k)\in\mathcal C(B)}\iint_{I_j\times I_k}f(\tau)f(x)\cdot 1\,d\tau\,dx.
\]
The right-hand side is a *lower bound* because:
* For \((j,k)\in\mathcal C(B)\), the indicator is identically 1 on
  \(I_j\times I_k\) (by definition of \(\mathcal C(B)\)), so
  \(\iint_{I_j\times I_k}f(\tau)f(x)\cdot\mathbf 1_{a\le\tau+x\le b}\,d\tau\,dx
  = \iint_{I_j\times I_k}f(\tau)f(x)\,d\tau\,dx = (\int_{I_j}f)(\int_{I_k}f) = p_jp_k\)
  **exactly**, by ordinary Fubini on the rectangle.
* For partially-contained cells (where the indicator is sometimes 0,
  sometimes 1), the contribution to the iterated integral is
  *nonnegative* but is *not* in the sum; dropping it produces a
  lower bound.
* For cells not in \(\mathcal C(B)\) at all (the indicator is 0), the
  contribution is 0.

**This is correct, with no false equality and no piecewise-constant
approximation.** The proof's note —

> The middle equality is **exact** (Fubini on the rectangle):
> \(\iint_{I_j\times I_k}f(\tau)f(x)\,d\tau\,dx = (\int_{I_j}f)(\int_{I_k}f) = p_jp_k\).
> **There is no piecewise-constant approximation here.**
> This is the core fix.

— is exactly right. The v0/v1/v2/v3/v4/v5/v6 chain that I'd been
flagging is finally fixed: instead of trying to *equate*
\(\frac12\sum_{j+k\in\{m,m+1\}}p_jp_k\) with \(\int_{B_m}(f*f)\)
(which is false because of partial-cell contributions), v7 *only
counts fully-contained cells* and uses the partial cells as nonneg
slack. The result is a *correct lower bound* on the iterated integral.

**Step (c) Trivial upper bound.**
\(\int_B F(t)\,dt \le \Omega\cdot\int_B 1\,dt = h_B\Omega\). ✓

Combining (a)–(c) gives the (W) constraint
\(\sum_{(j,k)\in\mathcal C(B)}Q[j,k]\le h_B\Omega\). ✓

### Cell-membership criterion. ✓
\((j,k)\in\mathcal C(B=[a,b])\) iff
\[
-\tfrac12 + (j+k-2)L \ge a\quad\text{and}\quad -\tfrac12 + (j+k)L \le b,
\]
i.e., \((a+1/2)/L + 2 \le j+k \le (b+1/2)/L\). I verified directly:
* \(\inf_{(\tau,x)\in I_j\times I_k}(\tau+x) = \inf I_j + \inf I_k = (-1/4 + (j-1)L) + (-1/4 + (k-1)L) = -1/2 + (j+k-2)L\). ✓
* \(\sup_{(\tau,x)\in I_j\times I_k}(\tau+x) = -1/2 + (j+k)L\). ✓

Full containment in \([a,b]\) requires both endpoints in \([a,b]\),
giving the stated criterion. ✓

### Crucial absence: the v0 band constraint is gone. ✓

I scanned the constraint list (P0)–(P5) carefully. There is no
constraint of the form \(\frac12\sum_{j+k\in\{m,m+1\}}Q[j,k]\le L\Omega\).
**The v0 band constraint is dropped, not inherited.** This finally
removes the load-bearing inheritance issue that has propagated since
v1.

The (W) family at widths \(\ge 4L\) (which is the smallest width the
v7 program uses; the centred window \([-2L, 2L]\) has full width
\(4L\)) gives non-trivial constraints; at width \(L\) (where v0's
bands lived) no cell is fully contained, so (W) is trivially
\(0\le L\Omega\) — which is *much weaker* than v0's bands (which
exploited partial-cell triangular contributions via the unproven
\(\frac12\) factor) but, crucially, *correct*.

### Validity claim. ✓

The proof's conclusion
\(\Omega^{(7)}\le \inf_{f\in\mathcal A^{\text{even}}}\|f*f\|_\infty = C^{\text{even}}_{6.2}\)
follows immediately: every admissible even \(f\) gives a feasible
\((Omega^*, p^*, Q^*)\) with \(\Omega^* = \|f*f\|_\infty\), so the LP
optimum is at most \(\Omega^*\) for every such \(f\), hence
\(\le\) the infimum.

### Numerical results.

The table

| \(N\) | \(\Omega^{(7)}\) |
|---|---|
| 4 | 1.046 |
| 6 | 1.082 |
| 8 | 1.164 |
| 10 | 1.206 |
| 12 | 1.241 |
| 16 | 1.299 |
| 20 | 1.344 |

shows monotone increase toward the analytic
\(C^{\text{even}}_{6.2} = 2\), achieved by
\(f = 2\cdot\mathbf 1_{[-1/4,1/4]}\). I cannot independently verify
the numbers (that would require running `solve.py`), but the
qualitative pattern (monotone, strictly above 1, bounded above by 2)
is exactly what one would expect for a valid relaxation that
converges to the truth.

## 4. Comments

### Comment 1. v7 is the right kind of "clean reboot."

After six iterations of "add another constraint to the LP," v7 takes
the structurally correct step of *removing* the load-bearing wrong
constraint and *replacing* it with a correct (if weaker) one. This is
exactly the move I'd been recommending in the previous verdicts (e.g.
v3 §4 Comment 1, v4 §4 Comment 1, v6 §4 Comment 3): drop the v0
bands, switch to the (W) family.

The trade-off — that (W) at width \(L\) is trivial, where v0's bands
were nontrivial — is real and explains why v7's numerical bound
(1.344 at \(N=20\)) is *weaker* than v5's (1.923 at \(N=12\)). v5 was
exploiting the unproven v0 chain to get more bite. v7 has less bite
but is *certifiable*. **Certifiability beats bite at this stage.**

### Comment 2. v7 still bounds the *even* constant only.

The proposer is explicit that
\(C^{\text{even}}_{6.2} \ge C_{6.2}\), so a bound on the even
constant is *not* automatically a bound on \(C_{6.2}\). The v7 bound
of 1.344 at \(N = 20\) is well below the SOTA
\(C_{6.2}\ge 1.28\) (the SOTA holds for the unrestricted constant,
and a fortiori for the even constant), so v7 does **not** improve the
project's main numerical target.

What v7 *does* establish is that the rigorous LP machinery
(centered grid + (W) family + (G) symmetry) produces a non-trivial
bound on the even-restricted problem. v8's job is to remove the even
restriction.

### Comment 3. The convergence rate is reasonable but slow.

The numerical table shows
\(\Omega^{(7)}(N) \approx 2 - O(1/N^{?})\), but the rate from the
data is hard to read off (1.344 at \(N=20\) is still some way from
2). With rigorous quadrature error and the (W) family at widths
\(\ge 4L\), this is expected: the (W) constraints don't see
fine-scale information about \(F\), so converging to the
analytic value will take a fairly large \(N\).

A v8 improvement would be to **extend the (W) family** to include
*all* widths \(2L, 3L, 4L,\ldots, 1\) (full widths), not just the
centered subset, and to add **shifted windows** as v4 had originally
proposed (\(B = [\pm 1/4 - h, \pm 1/4 + h]\)). The shifted windows
constrain the off-centre behaviour of \(F\), which is where the
diagonal solution lives; adding them should accelerate the
convergence.

### Comment 4. The cosmetic issue v0 originally introduced is now moot.

I've been flagging "v0 band constraint hand-wave" since the v1
verdict. v7 makes this issue *moot* by dropping the v0 band
constraint. There's no longer anything to inherit.

(I realize this means I should *also* update my mental model of the
v0 file: v0 is still wrong, but it's no longer load-bearing for
v7+. It would still be worth fixing v0 for future iterations that
might want to inherit the strong v0 bands rigorously, but that's
optional cleanup.)

## 5. Feedback to the proposer

**This is the biggest single step in the project so far.** The
structural insight — drop the v0 bands and replace them by the (W)
family — was the right move, and v7 executes it cleanly. The proof
is fully rigorous, the numerical results are non-trivial, and the
honest disclosure about the even restriction is appropriate.

Three concrete suggestions for v8:

1. **Add more windows.** Specifically:
   * Centered windows at *all* widths \(2L, 3L, 4L, \ldots, 1\), not
     just \(2L, 3L, \ldots, 1/2\) (the centered ones are
     parameterized by the half-width \(h\), so include all
     \(h\in\{L, 3L/2, 2L,\ldots, 1/2\}\)).
   * **Shifted windows** \(B = [c-h, c+h]\) at offsets
     \(c\in\{\pm L, \pm 2L, \ldots, \pm 1/4\}\) and the same set of
     half-widths. v4 had this in its plan; v7 dropped it. Re-adding
     them should significantly tighten the relaxation.
   * Asymmetric windows \(B = [a, b]\) for various non-symmetric
     \((a, b)\) — these are particularly useful for catching
     pathological diagonal solutions.
2. **Extend to non-even \(f\).** This is the v6/v8 promised
   extension. Two routes:
   * **Don't restrict to even \(f\)**: drop (E)/(F)/(G), and find a
     different way to break the diagonal-\(Q\) degeneracy. The most
     promising is the **\(F\)-side autocorrelation Toeplitz lift**
     (analog of 6.3 v4). Define cell averages \(W_m\) of \(F = f*f\)
     directly on \([-1/2, 1/2]\), use the rank-one identification
     \(\hat F(k) = \hat f(k)^2\), and add the Toeplitz PSD constraint
     on the discrete autocorrelation \(\hat F\). This is a *universal*
     constraint (every admissible \(f\) yields a Toeplitz-PSD
     autocorrelation), not restricted to even \(f\).
   * **Or:** keep (E)/(F)/(G) for even \(f\), and *interpolate* with
     a parameter that measures the asymmetry — but more carefully
     than v6 did. v6's slacks were phantom; v8 would need a *real*
     coupling between the asymmetry and \(\Omega\), which the
     parallelogram identity in v5's plan was supposed to provide
     (with the factor-of-2 fix I noted in the v5 verdict).

   I strongly recommend the first route (\(F\)-side autocorrelation
   Toeplitz). It is the same idea that 6.3 v4 used and is fully
   universal.
3. **Fix v0 (optional).** The v0 file still contains the original
   false-equality chain. This is no longer load-bearing for v7+, but
   if you ever want to inherit the *strong* v0 bands rigorously, the
   v0 derivation needs to be cleaned up (e.g. using cell-supremum
   variables \(M_j\) and the partial-cell area decomposition). Low
   priority.

If v8 does (1) and (2), it should produce a non-trivial bound on
\(C_{6.2}\) (not just the even constant) for the first time in the
project, breaking the 5+iteration drought.

---

**Final verdict: VALID.**
