# Verdict on v5/rigorousproof.md

## 1. Summary of the proof

v5 keeps everything from v4 and adds one more LMI on the same
$(R+1)\!\times\!(R+1)$ Toeplitz expression: the *dual* Bochner
constraint coming from $\Omega-M\ge 0$ on $[-2,2]$, namely
$$
\Omega I_{R+1}\;-\;T_{R+1}\;\succeq\;0. \tag{$\star\!\star\!\star\!\star\!\star$}
$$
Together with v4's $T_{R+1}\succeq 0$ this gives the spectral sandwich
$0\preceq T\preceq \Omega I$. The argument is the same trigonometric-
polynomial pairing as v4, applied to $\Omega-M$ in place of $M$;
the constant $\Omega$ contributes only on the diagonal because
$\int_{-2}^{2}e^{i\pi jx/2}dx = 4\delta_{j,0}$ for integer $j$.

## 2. Line-by-line verification

**Step 1 — $\Omega-M\ge 0$ on $[-2,2]$.**
By definition $\Omega:=\lVert M\rVert_\infty\ge M(x)$ a.e.\ on
$[-2,2]$. $M$ is continuous so the inequality holds pointwise on
$[-2,2]$. **No issue.**

**Step 2 — Bochner on $\Omega-M$.**

Independent algebraic check:
$\widehat{\Omega-M}(j)=\tfrac14\!\int_{-2}^{2}e^{-i\pi jx/2}(\Omega-M)dx$.
- $j=0$: $\tfrac14(4\Omega-1)=\Omega-\tfrac14$.
- $j\ne 0$: $\Omega\cdot\tfrac14\int e^{-i\pi jx/2}dx=0$ (integer $j$),
  leaving $-\hat M(j)$.

Therefore the Toeplitz matrix of $\Omega-M$ equals
$\Omega I_{R+1}-T_{R+1}$ entry-by-entry. Bochner on the
non-negative integrand $\Omega-M$ gives this matrix $\succeq 0$.
The proposer's algebra arrives at the same conclusion via the
quadratic-form route:
$$
0\;\le\;\int_{-2}^{2}p(x)(\Omega-M(x))\,dx
=4\!\sum_{k,l}u_k\overline{u_l}\bigl(\Omega\delta_{kl}-\overline{\hat M(k-l)}\bigr),
$$
where $p(x)=|\sum_k u_k e^{i\pi kx/2}|^{2}\ge 0$. The intermediate
identities
$\int_{-2}^{2}e^{i\pi(k-l)x/2}dx=4\delta_{kl}$ and
$\int_{-2}^{2}e^{i\pi jx/2}M\,dx=4\overline{\hat M(j)}$ are both
correct (the second was checked in v4). Conjugating the real sum
yields $u^{*}(\Omega I-T_{R+1})u\ge 0$ for all $u$, hence
$\Omega I-T_{R+1}\succeq 0$. **No issue.** The proposer's parenthetical
"wait, we only need…" interjection is conceptually right and not a
gap: Bochner needs $\Omega-M\ge 0$ only on the support of the
integral, and that holds on $[-2,2]$.

**Step 3 — validity at admissible $f$.**

Inherit the v4 admissible extension ($\tilde A_k=A_k^{\text{true}}$,
$b_k^{\text{prog}}=b_k^{\text{true}}$). Then $T_{R+1}$ in the program
becomes the true $[\hat M(k-l)]$, and Bochner on $\Omega-M\ge0$ gives
$\Omega I-T_{R+1}\succeq 0$. All other v0–v4 constraints continue to
hold as previously verified. So every admissible $f$ extends to a
v5-feasible point and the program optimum lower-bounds
$\inf_f\lVert M\rVert_\infty$. **No issue.**

**Step 4 — Convexity.**

$\Omega I_{R+1}-T_{R+1}$ has entries
- diagonal: $\Omega-\tfrac14$ — affine in $\Omega$;
- off-diagonal: $-(\tilde A_{|k-l|}-i\sigma(k-l)B_{|k-l|})/2$ —
  affine in $\tilde A,b^{\text{prog}}$.

Hermitian: $(\Omega I-T)_{lk}=\overline{(\Omega I-T)_{kl}}$
(inherited from $T$'s Hermitian-ness). So the constraint is a
proper LMI on an affine Hermitian matrix-valued map — convex. The
real-block reformulation is the same standard trick as v4 with the
same size. ✓

## 3. Strictness vs.\ v4

v5 = v4 + the new LMI $\Omega I-T\succeq 0$. Strictly more
constraints, so v5-feasible ⊂ v4-feasible. v5-optimum
$\Omega_{\text{v5}}^{\star}\ge\Omega_{\text{v4}}^{\star}$. The new
LMI is binding at any v4-optimum whose Toeplitz expression
$T^{(\text{v4-opt})}$ has $\lambda_{\max}(T^{(\text{v4-opt})}) > \Omega^{\text{v4-opt}}$
— generic unless v4 already happens to land in the bounded-density
cone.

## 4. Comments and feedback

* **Conceptually complete.** The pair $0\preceq T\preceq \Omega I$
  is the **standard SDP-expressible content** of the constraint
  "$M$ is a measure on $[-2,2]$ with bounded density in $[0,\Omega]$"
  truncated to the first $R+1$ Fourier modes. After v5 the
  Fourier-truncation Toeplitz programme is essentially saturated;
  further gains at this scale must come from larger $R$, longer
  cosine/sine tails (larger $T$), or *different* operators
  (e.g.\ moment matrices in $x$, or matrix Toeplitz blocks coupling
  $M$ with $M(\cdot+s)$ for translations $s$).

* **Specific corollaries.** The proposer's per-coefficient bound
  $A_k^{2}+B_k^{2}\le 4(\Omega-\tfrac14)^{2}$ from the $2\times 2$
  principal minors is a clean consequence. For the Erdős regime
  with $\Omega\approx 0.379$ this gives
  $A_k^{2}+B_k^{2}\lesssim 4\cdot(0.129)^{2}\approx 0.0666$,
  comparable to (but generically tighter than) what v1's energy
  bound gives when one $k$ dominates.

* **Implementation.** The dual LMI is the same data type as v4's
  PSD LMI, just shifted by $\Omega I$. From a CVXPY/MOSEK point of
  view the two LMIs can be combined into one block-diagonal LMI of
  size $2(R+1)$, halving the SDP setup overhead.

* **Caveat about $\Omega$ tightness.** The new constraint uses
  $\Omega$ as a free program variable, not as the *true*
  $\lVert M\rVert_\infty$. At the relaxation optimum the program's
  $\Omega$ is exactly the lower bound output; the LMI then says
  "the lifted Toeplitz expression must be compatible with a
  density bounded by that output $\Omega$." This is the right
  *direction* for a lower bound: a feasible $(c,d,\tilde A,B,\Omega)$
  with smaller $\Omega$ requires a Toeplitz expression squeezed
  into $[0,\Omega I]$, which is harder.

* **No new mathematical content beyond Bochner.** v5 is exactly the
  "second half" of the same theorem v4 applied — but applying
  Bochner to $\Omega-M$ rather than $M$ is genuinely a new
  *hypothesis* (the upper bound on $M$), as the proposer correctly
  highlights in §6. The proof is essentially identical to v4's,
  which is a feature, not a bug.

## 5. Final verdict

**VALID**
