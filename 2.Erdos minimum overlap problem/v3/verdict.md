# Verdict on v3/rigorousproof.md

## 1. Summary of the proof

v3 reintroduces the previously-discarded odd-index $A_{2r-1}^{2}$ family
into the Parseval‑energy bound by **lifting**: it adds a fresh program
variable $\tilde A_{2r-1}\in\mathbb R$ for each $r\in\{1,\dots,R\}$ and
constrains it with two convex one-sided inequalities:

* (A1) — Lemma 5 two-sided bound on the cosine Fourier coefficient via
  the $\alpha^\pm$ arrays;
* (A2) — Lemma 3 (3.5) upper bound:
  $\tilde A_{2r-1}\le\tfrac{4(-1)^{r+1}}{k\pi}a_k-2(a_k^{2}+b_k^{2})$
  (hypograph of a concave function).

The augmented energy constraint is

$$
\sum_{m=1}^{T}(c_m^{2}+d_m^{2})^{2}\;+\;4\sum_{r=1}^{R}\tilde A_{2r-1}^{2}\;+\;\frac{64}{\pi^{2}}\!\!\!\sum_{\substack{k=1\\k\text{ odd}}}^{2R-1}\!\!\!\frac{b_k^{2}}{k^{2}}\;\le\;2\Omega-\tfrac12. \tag{$\star\!\star\!\star$}
$$

Validity is shown by extending any admissible $f$ to a program point
with $\tilde A_{2r-1}^{\text{(adm)}}:=A_{2r-1}^{\text{true}}$, which
satisfies (A1) with margin, (A2) with equality, and reduces
$(\star\!\star\!\star)$ to the v2-style Parseval restriction with the
$A_{2r-1}^{2}$ block put back.

## 2. Line-by-line verification

**Step 0 (admissible extension).** Same scaffolding as v1/v2,
verified earlier. The new piece — $\tilde A_{2r-1}^{\text{(adm)}}:=A_{2r-1}$ —
is a definition; nothing to check yet. ✓

**Step 1 — (A1) at admissible $f$.**
White's Lemma 5 (lines 395–399 of `white_paper.txt`) is exactly the
two-sided bound used here:
$\tfrac{L}{2}\sum_{j}\alpha^{-}_{j,m}(w_j+v_j)\le A_m\le\tfrac{L}{2}\sum_{j}\alpha^{+}_{j,m}(w_j+v_j).$
The §5 choices
$\alpha^\pm_{j,m}=\cos(\pi mL(j-\tfrac12)/2)\pm\pi mL/4$ satisfy (3.9)
because by the mean-value theorem and
$|\tfrac{d}{dx}\cos(\pi mx/2)|=\tfrac{\pi m}{2}|\sin(\pi mx/2)|\le\tfrac{\pi m}{2}$,
the variation of $\cos(\pi mx/2)$ on an interval of length $L$ is at
most $\pi mL/4$. With $\tilde A_k=A_k$, (A1) is satisfied. **No issue.**

**Step 2 — (A2) at admissible $f$.**
White's (3.5) for odd $k=2r-1$: $\sin(k\pi/2)=(-1)^{r+1}$, so
$A_k=\tfrac{4(-1)^{r+1}}{k\pi}a_k-2(a_k^{2}+b_k^{2})$ holds as an
*equality*. Therefore (A2) is satisfied with slack 0. **No issue.**

**Step 3 — $(\star\!\star\!\star)$ at admissible $f$.**
From v2/v1: $\sum_{k\ge1}(A_k^{2}+B_k^{2})\le(\Omega-\tfrac14)/2$.
Discarded terms in (4) are
$\{A_{2m}^{2}+B_{2m}^{2}:m>T\}$ and $\{A_{2r-1}^{2}+B_{2r-1}^{2}:r>R\}$,
plus $B_{2m}^{2}\equiv 0$ — all $\ge 0$. The retained inequality is
$$
\sum_{m=1}^{T}A_{2m}^{2}+\sum_{r=1}^{R}(A_{2r-1}^{2}+B_{2r-1}^{2})\le\tfrac{\Omega-1/4}{2}.\tag{4}
$$
Substituting $A_{2m}^{2}=\tfrac14(c_m^{2}+d_m^{2})^{2}$,
$A_{2r-1}^{2}=(\tilde A_{2r-1}^{\text{(adm)}})^{2}$, and
$B_{2r-1}^{2}=\tfrac{16}{(2r-1)^{2}\pi^{2}}b_{2r-1}^{2}$, then
multiplying by 4 and re-indexing gives $(\star\!\star\!\star)$. **No issue.**

**Convexity.**
(A1): two affine inequalities in $(\tilde A_k,w,v)$ — convex. ✓
(A2): RHS is concave in $(a_k,b_k)$ (affine minus
$2\|(a_k,b_k)\|_2^{2}$); the inequality $\tilde A_k\le\text{concave}$
defines the hypograph of a concave function, which is convex. ✓
$(\star\!\star\!\star)$: $\|(c_m,d_m)\|_2^{4}$ convex (v1);
$\tilde A_{2r-1}^{2}$ univariate convex; $b_{2r-1}^{2}=(\text{affine expr})^{2}$
convex; sum with positive coefficients convex; RHS affine. ✓

**Lifting argument.**
The crucial conceptual point is that $\tilde A_k$ is *not* required to
equal $A_k^{\text{true}}$ at a generic feasible point — the program
merely keeps $\tilde A_k$ in the convex set
$\{(\tilde A_k,a_k,b_k,w,v):\text{(A1) and (A2)}\}$, which always
*contains* the admissible value $A_k^{\text{true}}$. So:

* **Admissibility ⇒ feasibility.** Every admissible $f$ extends to a
  program point with $\tilde A_k=A_k^{\text{true}}$, satisfying (A1),
  (A2), and $(\star\!\star\!\star)$.
* **Feasibility need not imply admissibility.** A generic feasible
  point picks $\tilde A_k$ inside its convex cell, possibly closer to
  $0$, which only *weakens* the energy constraint. The
  lower‑bound property is preserved because the program optimum
  $\Omega^{\star}\le \inf_f\lVert M\rVert_\infty$ requires *one*
  direction (every admissible $f$ is feasible), which is shown.

**No issue.**

## 3. Comments and feedback

* **Effectiveness.** This proposal is qualitatively the most ambitious
  step in the v0→v3 chain: it finally closes the energy decomposition
  with the *last* non-negative Parseval block. The trick — lifting the
  concave Lemma 3 RHS into a fresh DCP variable and recording only the
  *one-sided* inequalities that yield convex sets — is exactly the
  standard convex-hull / hypograph relaxation. Nicely executed.

* **Activeness.** As the proposal correctly notes, the new term
  $4\tilde A_{2r-1}^{2}$ is automatically inert when $0$ lies in the
  feasible cell (the program will set $\tilde A_{2r-1}=0$). The
  hard-square regime $p_1=p_2=0.381$, $-q_1=q_2=0.02$ forces
  $c_1\approx 0.381$ and $d_1\approx 0.02$. Then for $k=1$ the
  $\alpha^{-}$ side gives $F_{\cos}^{\text{LB}}[1]=A_2$-via-Lemma-5
  ≈ $-\tfrac12(c_1^{2}+d_1^{2})$ which is negative, while
  $F_{\cos}^{\text{UB}}[1]$ may be positive — so $0$ is generically
  inside the interval and the new term will *not* tighten the bound
  in this specific test case. **Caveat for the proposer:**
  the strictness gain may be hard to see in the validation regime;
  use a regime where one of the cosine‐sign brackets is forced
  positive (e.g.\ $p_1\to1/\pi$, $p_2\to1/\pi+\epsilon$) to observe
  the improvement empirically.

  Hmm — re-examining: at admissible $f$, $A_{2r-1}=A_{2r-1}^{\text{true}}$
  is *the* value of the cosine coefficient, not necessarily 0. The
  program nevertheless picks $\tilde A_{2r-1}$ to minimise its own
  contribution, which may differ from the admissible value. This
  divergence is innocuous for the *lower bound* validity, but it does
  mean the new constraint primarily binds when the **Lemma 5 envelope
  $[F^{\text{LB}},F^{\text{UB}}]$ is bounded away from 0** in a
  sign-definite way. In sign-indeterminate regimes the program just
  zeroes out the term.

* **Note on (A2) directionality.** (A2) records only the *upper*
  bound from (3.5). The *lower* bound on $\tilde A_{2r-1}$ from (A1)
  is via the $\alpha^{-}$ array; one could in principle add a second
  Lemma 3 inequality
  $\tilde A_{2r-1}\ge\tfrac{4(-1)^{r+1}}{k\pi}a_k-2(a_k^{2}+b_k^{2})$
  — but that would be the *lower* level set of a concave function,
  which is **non-convex** (the convex set is the hypograph, not the
  epigraph). The choice in v3 is the correct convex one. ✓

* **OCR caution.** As in v1/v2, the derivation re-uses (3.5) directly;
  I cross-checked (3.5) from Lemma 2 earlier.

## 4. Final verdict

**VALID**
