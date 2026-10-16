# A real proposer–verifier exchange (v9 → v10 → v14)

The early-version cluster of red ×'s in Figure 1 corresponds to a
multi-step exchange between the **proposer** (constructs SDP
relaxations and writes `rigorousproof.md`) and the **verifier**
(reads the proof and writes `verdict.md`). Below are excerpts from
three consecutive iterations that illustrate the loop.

The arc: v9 produces a strong-looking bound that only covers the
*even-restricted* constant $C^{\mathrm{even}}_{6.2}$; the verifier
accepts it under that caveat but flags the gap. v10 attempts a
non-even extension and collapses to $\Omega = 1$ due to a Fourier-side
*diagonal-inflation* degeneracy. The verifier diagnoses the
degeneracy and recommends a structural fix. Four iterations later,
v14 implements the fix and produces the **first rigorous lower bound
on the general $C_{6.2}$ in the project**.

---

## Panel 1 — Proposer (v9, `rigorousproof.md`)

> **v9 Theorem.**  $C^{\mathrm{even}}_{6.2} \ge \Omega^{(9)}$. At
> $(N, K) = (50, 64)$, $\Omega^{(9)} = 1.689$.
>
> The construction applies the Fejér kernel to $F = f * f$ for *even*
> $f$, using the identity
>
> $$F(t) = 1 + 2\sum_{k\ge 1} a_k^{2}\cos(2\pi k t), \qquad
>   a_k = \int_{-1/4}^{1/4} f(x)\cos(2\pi k x)\,dx,$$
>
> and the SDP lift $M[k,k] \ge a_k^{2}$ with $M \succeq 0$.

This is a strong bound, but **only for the even-restricted constant**.
The proposer is upfront about this: the next iteration will try to
drop the even restriction.

## Panel 2 — Verifier (`v10/verdict.md`)

> **VALID.**
>
> The math is correct: the constraints in v10 are mathematically
> valid (they hold for the rank-one lift from any admissible $f$),
> the Fejér convolution computation is right, the lifted constraints
> are well-defined.
>
> **However: the project does not advance numerically from v9 to v10.**
> v10's bound is $\Omega^{(10)} = 1$, strictly *worse* than v9's
> $\Omega^{(9)} = 1.689$ — because v10 doesn't have v9's
> even-symmetry constraints, and the Fourier-side machinery alone
> can't pin the degeneracy.
>
> **Diagnosis.** The PSD lift only forces *lower* bounds
> $M[k,k] \ge a_k^{2}$, $M[K+k, K+k] \ge b_k^{2}$. Nothing in the
> relaxation upper-bounds $M[K+k, K+k]$, so the solver inflates it
> freely, making the lifted constraint trivially satisfied.
>
> **Recommended fix.** Move to the *autocorrelation Toeplitz lift*:
> work with $|\hat F(k)|^{2}$ directly, where the manifest
> nonnegativity of $|\hat F(k)|^{2}$ blocks inflation.

The verifier accepts v10 as logically valid but flags it as
structurally vacuous, and prescribes a specific structural fix.

## Panel 3 — Proposer (v14, `rigorousproof.md`)

> **v14 Goal.**  Produce the *first rigorous lower bound on $C_{6.2}$
> in this project that does not require an even-$f$ restriction.*
>
> The construction follows the v9/v10 verdict's recommendation: apply
> the Fejér-kernel test functional to the **autocorrelation**
> $F * \check F$ (where $\check F(t) = F(-t)$) instead of $F$ itself.
>
> **Three properties of the autocorrelation:**
>
> (A) $(F * \check F)\widehat{}(k) = |\hat F(k)|^{2}$, real and
>     nonneg for every $k$.
>
> (B) $F * \check F \ge 0$ pointwise.
>
> (C) $\|F * \check F\|_\infty = (F * \check F)(0) = \int F^{2}\,dt$.
>
> **Key feature.** Every term in the lifted constraint has the
> same sign — there is no inflation degeneracy.
>
> **Result.** $\Omega^{(14)} = 1.063$ at $(N, K) = (16, 16)$.

This is the first teal point on the SOTA curve in Figure 1 — the
beginning of the climb that ends at $\Omega^{(51)} = 1.2937$.

---

## What this illustrates

* The verifier does **substantive mathematical work**, not just
  syntactic checking: in Panel 2 it correctly diagnoses *why* the
  SDP is vacuous (inflation of an un-upper-bounded diagonal) and
  proposes a *specific* structural remedy (autocorrelation Toeplitz).
* The proposer then **adopts the recommendation across versions** —
  v11–v13 are intermediate attempts at the fix, and v14 lands it.
* The arc is visible in Figure 1: v7–v13 are red ×'s (bound the wrong
  quantity, or vacuous), v14 is the first teal point on the running
  SOTA curve.
