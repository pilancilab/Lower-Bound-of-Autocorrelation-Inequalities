# v3 Proposal: convexify the missing odd-index $A_{k}^{2}$ family

## 1. Motivation

Constraint (v2) was
$$
\sum_{m=1}^{T}(c_m^{2}+d_m^{2})^{2}\;+\;\frac{64}{\pi^{2}}\sum_{\substack{k=1\\k\text{ odd}}}^{2R-1}\frac{b_k^{2}}{k^{2}}\;\le\;2\Omega-\tfrac12 \tag{v2}
$$
and was derived by retaining, in the energy Parseval bound
$\sum_{k\ge 1}(A_k^{2}+B_k^{2})\le(\Omega-1/4)/2$, only the two convex
families: even-index $A_{2m}^{2}$ and odd-index $B_{k}^{2}$. The one
remaining non-negative family is the odd-index $A_{2r-1}^{2}$. By White's
Lemma 3 (3.5),
$$
A_{2r-1}\;=\;\frac{4(-1)^{r+1}}{(2r-1)\pi}\,a_{2r-1}\;-\;2\bigl(a_{2r-1}^{2}+b_{2r-1}^{2}\bigr),
$$
a *concave* function of $(a_{2r-1},b_{2r-1})$; its square is generically
non-convex, which is why v2 dropped it.

The fix is a classical convexification. We introduce a *fresh program
variable* $\tilde A_{2r-1}\in\mathbb R$ and impose **only the directions of
the Lemma 3 / Lemma 5 inequalities that yield convex sets**. Then
$\tilde A_{2r-1}^{2}$ is the univariate square of a CVXPY variable —
convex — and can be added back into the Parseval bound.

White's program already uses Lemma 5's *lower* bound on $A_{m}$
(constraint (5.5)); the *upper* bound is silently dropped because it is
non-convex when $A_{m}$ is substituted by the concave Lemma 3
expression. Promoting $\tilde A_{2r-1}$ to a variable resolves this
exactly: every Lemma 5 / Lemma 3 inequality becomes a one-sided bound on
a scalar variable, all of them convex.

## 2. The new constraints

For every odd $k=2r-1\in\{1,3,\dots,2R-1\}$ introduce one fresh program
variable $\tilde A_{k}\in\mathbb R$. Impose

$$
\underbrace{\tfrac{L}{2}\textstyle\sum_{j=1}^{N}\alpha^{-}_{j,k}(w_{j}+v_{j})}_{=\,F_{\cos}^{\text{LB}}[k]}
\;\le\;\tilde A_{k}\;\le\;
\underbrace{\tfrac{L}{2}\textstyle\sum_{j=1}^{N}\alpha^{+}_{j,k}(w_{j}+v_{j})}_{=\,F_{\cos}^{\text{UB}}[k]}, \tag{A1}
$$
$$
\tilde A_{k}\;\le\;\frac{4(-1)^{r+1}}{k\pi}\,a_{k}\;-\;2\bigl(a_{k}^{2}+b_{k}^{2}\bigr). \tag{A2}
$$
Then strengthen the v2 energy inequality to

$$
\boxed{\;\sum_{m=1}^{T}(c_m^{2}+d_m^{2})^{2}\;+\;4\sum_{r=1}^{R}\tilde A_{2r-1}^{2}\;+\;\frac{64}{\pi^{2}}\sum_{\substack{k=1\\k\text{ odd}}}^{2R-1}\frac{b_k^{2}}{k^{2}}\;\le\;2\Omega-\tfrac12.\;} \tag{$\star\!\star\!\star$}
$$

## 3. Proof

We show that **(A1)+(A2)+($\star\!\star\!\star$) hold simultaneously at any
program point built from an admissible $(f,g,M)$ pair**. Combined with
the convexity argument in §4 this proves that the augmented convex
program remains a valid lower bound on $\inf_{f}\lVert M\rVert_{\infty}$.

### Step 0. Recall the admissible-to-program assignment.

Fix any admissible $f,g,M$ as in (2.1) of White. As in v1/v2, the
program-point assignment uses
$w_{j}=\tfrac1L\!\int_{(j-1)L}^{jL}\!M$,
$v_{j}=\tfrac1L\!\int_{-jL}^{-(j-1)L}\!M$,
$c_{k},d_{k}$ are the (3.1) Fourier coefficients of $f$ on $[-1,1]$,
$a_{k},b_{k}$ are the (3.2) Fourier coefficients of $f$ on $[-2,2]$,
$\varepsilon_{k},\delta_{k}$ capture the Lemma 4 tails so that the
program's $a_{k}^{\text{prog}},b_{k}^{\text{prog}}$ equal the true
$a_{k},b_{k}$ of $f$.

We **additionally** set
$$
\tilde A_{k}^{\text{(adm)}}\;:=\;A_{k}\;=\;\frac{4\sin(k\pi/2)}{k\pi}\,a_{k}-2\bigl(a_{k}^{2}+b_{k}^{2}\bigr)\quad(k=2r-1\in\{1,3,\dots,2R-1\}),
$$
i.e.\ the true value of the cosine Fourier coefficient of $M$ on
$[-2,2]$.

### Step 1. (A1) holds.

By Lemma 5 of White (combined with our choice
$\alpha^{-}_{j,k}=\cos(\pi k L(j-\tfrac12)/2)-\tfrac{\pi kL}{4}$ and
$\alpha^{+}_{j,k}=\cos(\pi k L(j-\tfrac12)/2)+\tfrac{\pi kL}{4}$, both
satisfying the inequalities (3.9) because
$\lvert\!\tfrac{d}{dx}\cos(\pi kx/2)\rvert\le \tfrac{\pi k}{2}$),
$$
\tfrac{L}{2}\sum_{j=1}^{N}\alpha^{-}_{j,k}(w_{j}+v_{j})\;\le\;A_{k}\;\le\;\tfrac{L}{2}\sum_{j=1}^{N}\alpha^{+}_{j,k}(w_{j}+v_{j}).
$$
With $\tilde A_{k}^{\text{(adm)}}=A_{k}$, this is (A1). ✓

### Step 2. (A2) holds.

For $k=2r-1$, $\sin(k\pi/2)=(-1)^{r+1}$, so (3.5) gives
$A_{k}=\tfrac{4(-1)^{r+1}}{k\pi}a_{k}-2(a_{k}^{2}+b_{k}^{2})$, i.e.\
*equality* in (A2). Hence (A2) is satisfied with slack 0. ✓

### Step 3. ($\star\!\star\!\star$) holds.

From v1/v2 we already have the global Parseval bound
$$
\sum_{k=1}^{\infty}(A_{k}^{2}+B_{k}^{2})\;\le\;\tfrac{\Omega-1/4}{2}. \tag{3}
$$
Discard only the strictly non-negative tail terms outside the program's
truncation horizons; every dropped term is $\ge0$:

* $A_{2m}^{2}+B_{2m}^{2}=A_{2m}^{2}\ge0$ for $m>T$ (no $c_m$ in program);
* $A_{2r-1}^{2}+B_{2r-1}^{2}\ge0$ for $r>R$ (no $b_{2r-1}$ in program);
* $B_{2m}^{2}=0$ (Lemma 3 (3.6) with $\sin(m\pi)=0$).

This leaves
$$
\sum_{m=1}^{T}A_{2m}^{2}\;+\;\sum_{r=1}^{R}\bigl(A_{2r-1}^{2}+B_{2r-1}^{2}\bigr)\;\le\;\tfrac{\Omega-1/4}{2}. \tag{4}
$$

Substituting:

* $A_{2m}^{2}=\tfrac14(c_m^{2}+d_m^{2})^{2}$ (v1, Lemma 3 even branch);
* $A_{2r-1}^{2}=(\tilde A_{2r-1}^{\text{(adm)}})^{2}$ by definition of
  $\tilde A_{2r-1}^{\text{(adm)}}$ in Step 0;
* $B_{2r-1}^{2}=\tfrac{16}{(2r-1)^{2}\pi^{2}}b_{2r-1}^{2}$ (v2, Lemma 3 (3.6)).

(4) becomes
$$
\tfrac14\sum_{m=1}^{T}(c_m^{2}+d_m^{2})^{2}\;+\;\sum_{r=1}^{R}\tilde A_{2r-1}^{2}\;+\;\frac{16}{\pi^{2}}\sum_{r=1}^{R}\frac{b_{2r-1}^{2}}{(2r-1)^{2}}\;\le\;\frac{\Omega-1/4}{2}.
$$
Multiplying by 4 and re-indexing $k=2r-1$ for the odd-$b$ sum
gives $(\star\!\star\!\star)$. ✓

Therefore every admissible $f$ extends to a feasible point of the
augmented program (with $\tilde A_{k}$ set to the true $A_{k}$), so the
program's optimum continues to be a valid lower bound on
$\inf_{f}\lVert M\rVert_{\infty}$. $\blacksquare$

## 4. Convexity of all new constraints

* **(A1).** $F_{\cos}^{\text{LB}}[k]$ and $F_{\cos}^{\text{UB}}[k]$ are
  *linear* functions of $(w_j,v_j)$ (the $\alpha^{\pm}_{j,k}$ are
  constants). Hence (A1) is two affine inequalities in
  $(\tilde A_k,w,v)$ — convex.

* **(A2).** The right-hand side
  $\tfrac{4(-1)^{r+1}}{k\pi}a_{k}-2(a_{k}^{2}+b_{k}^{2})$ is the sum of
  a linear (in $a_k$) and a concave (negative norm-squared in
  $(a_k,b_k)$) term — globally **concave** in $(a_k,b_k)$. The
  inequality $\tilde A_k\le\text{concave}(a_k,b_k)$ defines the
  *hypograph* of a concave function, which is a convex set in
  $(\tilde A_k,a_k,b_k)$. ✓

* **($\star\!\star\!\star$).**
  $(c_m^{2}+d_m^{2})^{2}=\lVert(c_m,d_m)\rVert_{2}^{4}$ is convex (v1);
  $\tilde A_{2r-1}^{2}$ is a univariate square of a CVXPY variable —
  convex; $b_{2r-1}^{2}$ is the square of an affine program expression
  (the truncated Lemma 4 representation of $b_{2r-1}$) — convex.
  All multipliers in the LHS are positive constants; sum of convex is
  convex. RHS $2\Omega-\tfrac12$ is affine in $\Omega$. ✓

Implementation-wise, the v2 slacks $s_m\ge c_m^{2}+d_m^{2}$ are reused;
each $\tilde A_{2r-1}^{2}$ becomes either a `cp.square` of a scalar or
an extra rotated SOC; the whole augmented LHS is a sum of squares of
affine quantities — a single second-order cone constraint with the
same shape as v2 but with $R$ extra terms.

## 5. Why this is genuinely new

* The Lemma 5 *upper* bound on the cosine Fourier coefficients of $M$
  is, **for the first time**, used. (White's (5.5) only encodes the
  lower bound; the upper bound is non-convex when written directly
  through the concave Lemma 3 expression — but is perfectly convex with
  $\tilde A_k$ as a free variable.)
* The odd-index $A_{2r-1}^{2}$ family — the last non-negative term in
  the Parseval bound — is now included in a *convex* way, completing
  the energy-decomposition programme begun in v1.
* Strictly stronger than v2: any feasible v3 point is feasible for v2
  (set $\tilde A_{2r-1}=0$ if $0$ lies in the interval
  $[F_{\cos}^{\text{LB}}[k],\,\min(F_{\cos}^{\text{UB}}[k],
  \text{Lemma 3 RHS})]$, else the closer endpoint), so v3's optimum
  $\Omega^{\star}\ge$ v2's optimum.

## 6. Note on the relaxation

The program does **not** force $\tilde A_{2r-1}$ to equal the true
$A_{2r-1}$; it only constrains it via (A1) and (A2). Hence the
program may choose $\tilde A_{2r-1}$ to *minimise* its contribution to
the LHS of ($\star\!\star\!\star$):

* If $0\in[F_{\cos}^{\text{LB}}[k],\,\min(F_{\cos}^{\text{UB}}[k],\,\text{Lemma 3 RHS})]$,
  the program sets $\tilde A_{2r-1}=0$ and the new term contributes
  nothing — v3 then coincides with v2.
* If $F_{\cos}^{\text{LB}}[k]>0$ (a sign-definite cosine), the program
  is forced to take $\tilde A_{2r-1}\ge F_{\cos}^{\text{LB}}[k]$ and
  $\tilde A_{2r-1}^{2}\ge(F_{\cos}^{\text{LB}}[k])^{2}>0$, adding a
  strictly positive contribution.

So the augmentation is **automatically inert in sign-indeterminate
regimes and automatically active in sign-determinate ones**, which is
the precise behaviour the verdict on v2 called for.
