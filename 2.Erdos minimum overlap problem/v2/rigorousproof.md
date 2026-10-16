# v2 Proposal: strengthen v1 with the odd-index $B_k^{2}$ terms

## 1. Motivation

The v1 constraint
$$
\sum_{m=1}^{T}(c_{m}^{2}+d_{m}^{2})^{2}\;\le\;2\Omega-\tfrac12 \tag{v1}
$$
is the inequality obtained from $\int_{-2}^{2}M^{2}\le \Omega$ after
restricting the Parseval sum $\sum_{k\ge1}(A_{k}^{2}+B_{k}^{2})$ to **only
the even-index $A_{k}^{2}$ terms**. Two non‑negative families were thrown
away:

* odd-index $A_{k}^{2}$ (these are quartic in $a_{k},b_{k}$ — **not convex**, so we still drop them);
* odd-index $B_{k}^{2}$ — but by Lemma 3 (3.6) of White,
  $B_{k}=-\tfrac{4\sin(k\pi/2)}{k\pi}\,b_{k}$, so for odd $k$
  $$
  B_{k}^{2}=\frac{16}{k^{2}\pi^{2}}\,b_{k}^{2},
  $$
  which is **convex** (a scaled square of the program variable $b_{k}$).

Therefore we may *keep* the odd-$k$ sine energy and obtain a strictly
stronger convex constraint than $(\text{v1})$.

## 2. The new constraint

**Theorem 2 (Strengthened quartic energy constraint).** *Let
$f,g,M$ satisfy White's (2.1), set $\Omega=\lVert M\rVert_\infty$, let
$\{c_k,d_k\}_{k\ge1}$ be the sine–cosine Fourier coefficients of $f$ on
$[-1,1]$, and let $\{b_{k}\}_{k\ge1}$ be the sine Fourier coefficients of
$f$ on $[-2,2]$ as in (3.2) of White. Then for all integers $T,R\ge1$:*
$$
\boxed{\;\sum_{m=1}^{T}(c_{m}^{2}+d_{m}^{2})^{2}\;+\;\frac{64}{\pi^{2}}\sum_{\substack{k=1\\k\text{ odd}}}^{2R-1}\frac{b_{k}^{2}}{k^{2}}\;\le\;2\Omega-\tfrac12.\;} \tag{$\star\star$}
$$

## 3. Proof

We reuse Steps 1–4 of the v1 proof verbatim (only Steps 5 and 6 change).

### Steps 1–4 (recap).

The pointwise bound $0\le M\le\Omega$ combined with $\int_{-2}^{2}M=1$
and Parseval on $[-2,2]$ gives
$$
\sum_{k=1}^{\infty}(A_{k}^{2}+B_{k}^{2})\;\le\;\tfrac{\Omega-1/4}{2}. \tag{3}
$$

### Step 5′. Discard only the **non-convex** non-negative terms.

In (3), every $A_{k}^{2}+B_{k}^{2}\ge 0$. We discard exactly the terms
that are *not convex in our program variables*:

* $A_{k}^{2}$ for **odd** $k$ — these are quartic in $(a_{k},b_{k})$
  with the wrong sign of the leading quadratic, so they are
  non-convex; drop them.
* $A_{k}^{2}+B_{k}^{2}$ for $k=2m$ with $m>T$ — the program has no
  variable $c_{m}$ in that range; drop them.
* $A_{k}^{2}+B_{k}^{2}$ for $k=2r-1$ with $2r-1>2R-1$, i.e.\ $r>R$ —
  the program has no variable $b_{k}$ in that range; drop them.
* $B_{k}^{2}$ for **even** $k$ — these are identically zero (Lemma 3,
  (3.6)), so dropping them costs nothing.

Every dropped term is non-negative, so:
$$
\underbrace{\sum_{m=1}^{T} A_{2m}^{2}}_{\text{even-index cosine energy}}
\;+\;\underbrace{\sum_{r=1}^{R} B_{2r-1}^{2}}_{\text{odd-index sine energy}}
\;\le\;\sum_{k=1}^{\infty}(A_{k}^{2}+B_{k}^{2})
\;\le\;\tfrac{\Omega-1/4}{2}. \tag{4′}
$$

### Step 6′. Evaluate the two surviving families via Lemma 3.

**Even cosine, $k=2m$:** exactly as in v1, equation (3.5) gives
$A_{2m}=-\tfrac12(c_{m}^{2}+d_{m}^{2})$, hence
$A_{2m}^{2}=\tfrac14(c_{m}^{2}+d_{m}^{2})^{2}$.

**Odd sine, $k=2r-1$:** by (3.6) of Lemma 3 in White,
$$
B_{2r-1}\;=\;-\frac{4\sin((2r-1)\pi/2)}{(2r-1)\pi}\,b_{2r-1}
\;=\;-\frac{4(-1)^{r+1}}{(2r-1)\pi}\,b_{2r-1}.
$$
The sign is irrelevant after squaring, and $\sin^{2}((2r-1)\pi/2)=1$, so
$$
B_{2r-1}^{2}\;=\;\frac{16}{(2r-1)^{2}\pi^{2}}\,b_{2r-1}^{2}. \tag{6}
$$

Substituting both into (4′):
$$
\tfrac14\sum_{m=1}^{T}(c_{m}^{2}+d_{m}^{2})^{2}
\;+\;\frac{16}{\pi^{2}}\sum_{r=1}^{R}\frac{b_{2r-1}^{2}}{(2r-1)^{2}}
\;\le\;\frac{\Omega-1/4}{2}.
$$

Multiplying through by 4 and re-indexing $k=2r-1$ for odd
$k\in\{1,3,\dots,2R-1\}$ gives exactly $(\star\star)$. $\blacksquare$

## 4. Strictness

Constraint $(\star\star)$ is strictly stronger than (v1) whenever any of
the program variables $b_{1},b_{3},\dots,b_{2R-1}$ is non-zero:
the new term $\tfrac{64}{\pi^{2}}\sum b_{k}^{2}/k^{2}$ is non-negative
and adds to the LHS. For $b_{k}=0$ on odd $k$ the two constraints
coincide.

In White's "hard square" regime
$h_{1}=h_{2}=0.015$, $p_{1}=p_{2}=0.381$, $-q_{1}=q_{2}=0.02$, $d_{1}$ is
small but non-zero, and the relation
$b_{1}=\delta_{1}+\tfrac{4}{\pi}\sum_{k=1}^{T}\tfrac{k(-1)^{k}\sin(\pi/2)}{1-4k^{2}}d_{k}$
ties $b_{1}$ to *all* $d_{k}$, so the new term is generically active.

## 5. Convexity

Each $b_{k}$ for odd $k\le 2R-1$ is a CVXPY scalar variable. $b_{k}^{2}$
is convex, the multiplier $\frac{64}{k^{2}\pi^{2}}>0$ is constant, and
$\sum$ of convex is convex. Combined with the v1 quartic term (already
shown convex), the LHS of $(\star\star)$ is convex; the RHS is affine in
$\Omega$. Hence $(\star\star)$ defines a convex feasible region. The
implementation reuses the v1 slacks $s_{m}\ge c_{m}^{2}+d_{m}^{2}$ for
the quartic part and simply adds the linear-quadratic odd-$k$ piece into
the same SOC.

## 6. Validity within the program: a careful note on $b_{k}$

The program does not have access to the "true" $b_{k}$; it stores the
**truncated** version
$$
b_{k}^{\text{prog}}\;=\;\delta_{k}\;+\;\frac{4}{\pi}\sum_{j=1}^{T}\frac{j(-1)^{j}\sin(k\pi/2)}{k^{2}-4j^{2}}\,d_{j},
$$
where $\delta_{k}$ is a free variable in $[-\Delta_{k},\Delta_{k}]$ with
$\Delta_{k}$ exactly the Lemma 4 tail bound. As White shows in the proof
of Proposition 9 of his paper, the assignment $\delta_{k}=$ (full tail
sum) makes $b_{k}^{\text{prog}}$ *equal* to the true Fourier coefficient
$b_{k}$ of any admissible $f$, and Lemma 4 certifies that this
assignment is feasible. Therefore *every* admissible $f$ yields a
feasible point of the program with $b_{k}^{\text{prog}}=b_{k}^{\text{true}}$,
and at such points $(\star\star)$ becomes the inequality proved in
Sections 2–3 above. Hence $(\star\star)$ — imposed as a *program*
constraint on the variable $b_{k}^{\text{prog}}$ — preserves all
admissible feasible points and the convex optimum remains a valid lower
bound on $\inf_{f}\sup_{x}M$.
