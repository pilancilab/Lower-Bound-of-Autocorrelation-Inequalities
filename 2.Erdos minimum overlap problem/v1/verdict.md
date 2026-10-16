# Verdict on v1/rigorousproof.md

## 1. Summary of the proof

The proposal introduces one new convex constraint linking the
Fourier‑energy of $f$ on $[-1,1]$ to $\Omega=\lVert M\rVert_\infty$:

$$
\sum_{m=1}^{T}\bigl(c_m^{2}+d_m^{2}\bigr)^{2}\;\le\;2\Omega-\tfrac12.
\tag{$\star$}
$$

The derivation is a six‑step chain:

1. Pointwise: $0\le M(x)\le \Omega$ implies $M(x)^{2}\le \Omega\,M(x)$.
2. Integrate over $[-2,2]$ and use $\int_{-2}^{2}M=1$ (White (2.2)),
   yielding $\int_{-2}^{2}M^{2}\le \Omega$.
3. Parseval on $[-2,2]$ with $M=\tfrac14+\sum(A_k\cos+B_k\sin)$ gives
   $\int_{-2}^{2}M^{2}=\tfrac14+2\sum_{k\ge1}(A_k^{2}+B_k^{2})$.
4. Combine: $\sum_{k\ge1}(A_k^{2}+B_k^{2})\le(\Omega-\tfrac14)/2$.
5. Discard all non‑$A_{2m}^{2}$ terms (each is $\ge 0$):
   $\sum_{m=1}^{T}A_{2m}^{2}\le(\Omega-\tfrac14)/2$.
6. Use White's Lemma 3, even branch: $a_{2m}=c_m/2$, $b_{2m}=d_m/2$ and
   $\sin(m\pi)=0$ kill the linear part of (3.5), giving
   $A_{2m}=-\tfrac12(c_m^{2}+d_m^{2})$, hence
   $A_{2m}^{2}=\tfrac14(c_m^{2}+d_m^{2})^{2}$.
   Multiplying through yields $(\star)$.

A short convexity note (sum of $\lVert(c_m,d_m)\rVert_2^{4}$ functions vs.
an affine RHS) and an SOC implementation outline complete the file.

## 2. Line-by-line verification

**Step 1 (pointwise quadratic bound).**
$M(x)=\int_{-1}^{1}f(t)g(x+t)\,dt$ with $f,g\in[0,1]$ implies $M\ge 0$.
$M$ is continuous (Lipschitz, in fact), so $\Omega:=\lVert M\rVert_\infty
=\sup_{x\in[-2,2]}M(x)$ is achieved and $M(x)\le\Omega$ pointwise.
The implication $0\le M\le\Omega\Rightarrow M(\Omega-M)\ge 0
\Leftrightarrow M^{2}\le\Omega M$ is trivial algebra on a real‑valued
function. **No issue.**

**Step 2 (integration).**
The proof imports "$f,g,M$ satisfy (2.1) of White". White's (2.1)
includes the normalisation $\int_{-1}^{1}f\,dx=1$, which is exactly
what makes White's (2.2) $\int_{-2}^{2}M\,dx=1$ valid. Under the same
hypothesis this step is legal. **No issue.**

(Caveat for the proposer: Problem.md does not literally state
$\int f=1$; it inherits it implicitly because the cited SOTA 0.379005
is White's bound. The proof keeps this implicit by reusing White's
(2.1) — that is fine, but the proposer should not later attempt to
weaken (2.1) without retracing this step.)

**Step 3 (Parseval on $[-2,2]$).**
With $A_0=\tfrac14$ (forced by $\int_{-2}^{2}M=1$) and the
orthogonality numbers $\int_{-2}^{2}1\,dx=4,\
\int_{-2}^{2}\cos^{2}(k\pi x/2)\,dx=\int_{-2}^{2}\sin^{2}(k\pi x/2)\,dx=2$,
$$
\int_{-2}^{2}M^{2}\,dx
   = (\tfrac14)^{2}\cdot 4 \;+\; 2\sum_{k\ge1}(A_k^{2}+B_k^{2})
   = \tfrac14 \;+\; 2\sum_{k\ge1}(A_k^{2}+B_k^{2}).
$$
$M\in L^{2}([-2,2])$ since $M$ is bounded, so Parseval applies.
**No issue.**

**Step 4 (combine).** Algebraic rearrangement. **No issue.**

**Step 5 (discard non‑negative tail).**
Every $A_k^{2}+B_k^{2}\ge 0$, so restricting the sum to $\{A_{2m}^{2}\}_{1\le m\le T}$
weakens it; the inequality is preserved. **No issue.**

**Step 6 (Lemma 3, even branch).**
I re‑derived White's (3.5) directly from Lemma 2 to neutralise the
OCR artefact in `white_paper.txt` (which renders the squared terms as
"a2m+b2m"). Using $A_m=\hat M(m)+\hat M(-m)$,
$2\hat f(\pm m)=a_m\mp i b_m$, $|\hat f(\pm m)|^{2}=(a_m^{2}+b_m^{2})/4$,
$$
A_m=\frac{4\sin(m\pi/2)}{m\pi}\,a_m\;-\;2(a_m^{2}+b_m^{2}).
$$
Substituting $m\mapsto 2m$ gives
$A_{2m}=\frac{4\sin(m\pi)}{2m\pi}a_{2m}-2(a_{2m}^{2}+b_{2m}^{2})
=-2(a_{2m}^{2}+b_{2m}^{2})$. ✓

From White's (3.3) and (3.4), even branch: $a_{2m}=\tfrac12 c_m$,
$b_{2m}=\tfrac12 d_m$. ✓

Therefore $A_{2m}=-\tfrac12(c_m^{2}+d_m^{2})$ and
$A_{2m}^{2}=\tfrac14(c_m^{2}+d_m^{2})^{2}$. ✓

Substituting into Step 5 and clearing the $\tfrac14$ gives $(\star)$.
**No issue.**

**Convexity (Section 4).**
$(c_m,d_m)\mapsto(c_m^{2}+d_m^{2})^{2}=\lVert(c_m,d_m)\rVert_2^{4}$.
$\lVert\cdot\rVert_2^{2}$ is convex on $\mathbb R^{2}$; $t\mapsto t^{2}$
is convex and non‑decreasing on $t\ge 0$; composition of a
non‑decreasing convex function with a convex function is convex
(Boyd–Vandenberghe, §3.2). Sums preserve convexity and the RHS
$2\Omega-\tfrac12$ is affine in $\Omega$. ✓

**SOC implementation (Section 6).**
With slacks $s_m\ge c_m^{2}+d_m^{2}$ (rotated SOC, DCP‑compliant) and
$\sum_{m}s_m^{2}\le 2\Omega-\tfrac12$ (a single SOC),
every admissible $f$ extends to a feasible tuple by choosing
$s_m=c_m^{2}+d_m^{2}$; hence the SOCP optimum is still a valid lower
bound for $\inf_f\sup_xM$. (Note that the slack relaxation is in fact
slightly tighter than $(\star)$, not looser, which is acceptable for a
lower bound because admissibility is preserved.) ✓

## 3. Comments and feedback

* **Novelty.** The constraint is genuinely new relative to White's
  program: it bounds the $\ell^{2}$‑norm of the sequence
  $(c_m^{2}+d_m^{2})_m$, whereas White's (5.11) is the $\ell^{1}$ bound
  $\sum(c_m^{2}+d_m^{2})\le\tfrac12$. Concretely, the test case
  $c_1^{2}+d_1^{2}=\tfrac12$ with all higher coefficients zero is
  feasible under (5.11) at $\Omega=\tfrac14$, but $(\star)$ forces
  $\Omega\ge\tfrac38=0.375$. So the new constraint actually has bite
  in regimes the existing program misses.

* **Sharpness.** The pointwise bound $M(\Omega-M)\ge 0$ is the
  *coarsest* non‑trivial pointwise relation. A natural next step is
  to replace it by $\int(M-\Omega)_+\,M=0$, or to use
  $(M-\alpha)(\Omega-M)\ge 0$ on the level set $\{M\ge\alpha\}$, which
  would give a parameterised family of quartic constraints and
  potentially a stronger Fourier bound on $\sum(A_k^{2}+B_k^{2})$.

* **Care with normalisations.** Future proofs should explicitly cite
  $\int_{-1}^{1}f=1$ wherever (2.2) is used, since Problem.md alone
  does not impose it.

* **Care with the white_paper.txt OCR.** The text extract renders
  $a_m^{2}+b_m^{2}$ as "a2m+b2m". Future proofs should re‑derive
  Lemma 3 from Lemma 2 (or from the source PDF) rather than copying
  the text file, to avoid propagating that artefact.

## 4. Final verdict

**VALID**
