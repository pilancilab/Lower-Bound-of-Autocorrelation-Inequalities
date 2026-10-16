# v1 Proposal: a new convex constraint from the $L^2$ energy of $M$

## 1. Motivation

White's program in §5 of *Erdős' minimum overlap problem* (arXiv:2201.05704)
encodes the following pieces of information about an admissible pair
$f,g$ with $g=1-f$ on $[-1,1]$ and $M(x)=\int_{-1}^{1}f(t)g(x+t)\,dt$:

* $0 \le M(x)\le \Omega$ pointwise (`(5.1)`);
* $\int_{-2}^{2} M = 1$ (`(5.2)`);
* first–moment / second–moment bounds (Lemma 7, `(5.3)–(5.4)`);
* Fourier cosine / sine identities and inequalities from Lemmas 3 and 5
  (`(5.5)–(5.7)`);
* Parseval bound $\sum_{k\ge1}(c_k^{2}+d_k^{2})\le 1/2$ (`(5.11)`).

The Parseval bound is just the **$L^{2}$ energy of $f$** translated to its
Fourier coefficients. The energy of $M$ on $[-2,2]$ — which is also an
infinite series in the Fourier coefficients $A_{k},B_{k}$ of $M$ — is
**never used**. We exploit it.

The new ingredient is the trivial pointwise inequality
$M(x)\bigl(\Omega - M(x)\bigr)\ge 0$, which gives
$\int_{-2}^{2} M^{2}\le \Omega \int_{-2}^{2}M = \Omega$. Translating this
through Parseval and Lemma 3 produces a convex quartic constraint that
ties the Fourier energy of $f$ directly to the bound $\Omega$ — something
no existing constraint of the program does.

## 2. The new constraint

**Theorem 1 (Quartic energy constraint).** *Let $f,g,M$ satisfy (2.1) of
White's paper, set $\Omega=\lVert M\rVert_\infty$, and let
$\{c_k,d_k\}_{k\ge1}$ be the sine–cosine Fourier coefficients of $f$ on
$[-1,1]$ as in (3.1). Then for every integer $T\ge 1$:*
$$
\boxed{\;\sum_{m=1}^{T}\bigl(c_{m}^{2}+d_{m}^{2}\bigr)^{2}\;\le\;2\,\Omega-\tfrac{1}{2}.\;}
\tag{$\star$}
$$

## 3. Proof

The proof has six self-contained steps. We use **only** White's notation
and his Lemmas 2 and 3.

### Step 1. Pointwise quadratic bound on $M$.

For every $x\in\mathbb{R}$ we have $0\le M(x)\le \Omega$ (with $M\equiv 0$
outside $[-2,2]$). Therefore
$$
M(x)\bigl(\Omega-M(x)\bigr)\;\ge\;0
\quad\Longleftrightarrow\quad
M(x)^{2}\;\le\;\Omega\,M(x).
$$

### Step 2. Integrate.

Integrating Step 1 over $[-2,2]$ and using
$\int_{-2}^{2}M(x)\,dx=1$ (equation (2.2) of White),
$$
\int_{-2}^{2} M(x)^{2}\,dx\;\le\;\Omega\!\int_{-2}^{2}M(x)\,dx\;=\;\Omega. \tag{1}
$$

### Step 3. Parseval on $[-2,2]$.

Recall the Fourier series (3.2) of White, written on the interval
$[-2,2]$ (period 4):
$$
M(x)=\tfrac14+\sum_{k=1}^{\infty}\bigl(A_{k}\cos(k\pi x/2)+B_{k}\sin(k\pi x/2)\bigr).
$$
The basis functions $1,\cos(k\pi x/2),\sin(k\pi x/2)$ for $k\ge1$ are
orthogonal on $[-2,2]$ with
$$
\int_{-2}^{2}1\,dx=4,\quad \int_{-2}^{2}\cos^{2}(k\pi x/2)\,dx=2,\quad \int_{-2}^{2}\sin^{2}(k\pi x/2)\,dx=2.
$$
Hence Parseval's identity gives
$$
\int_{-2}^{2} M(x)^{2}\,dx \;=\;\bigl(\tfrac14\bigr)^{2}\!\cdot 4+2\sum_{k=1}^{\infty}\bigl(A_{k}^{2}+B_{k}^{2}\bigr)\;=\;\tfrac14+2\sum_{k=1}^{\infty}\bigl(A_{k}^{2}+B_{k}^{2}\bigr). \tag{2}
$$

### Step 4. Combine (1) and (2).

Substituting (2) into (1):
$$
\tfrac14+2\sum_{k=1}^{\infty}(A_{k}^{2}+B_{k}^{2})\;\le\;\Omega
\;\Longleftrightarrow\;
\sum_{k=1}^{\infty}(A_{k}^{2}+B_{k}^{2})\;\le\;\tfrac{\Omega-1/4}{2}. \tag{3}
$$

### Step 5. Discard non‑negative tail terms.

Every term $A_{k}^{2}+B_{k}^{2}$ in (3) is non-negative. In particular,
restricting to the *even-indexed cosine* terms with $k=2m$ for
$1\le m\le T$,
$$
\sum_{m=1}^{T} A_{2m}^{2}\;\le\;\sum_{k=1}^{\infty}(A_{k}^{2}+B_{k}^{2})\;\le\;\tfrac{\Omega-1/4}{2}. \tag{4}
$$

### Step 6. Evaluate $A_{2m}$ via Lemma 3.

By formula (3.5) of Lemma 3,
$$
A_{2m}=\frac{4\sin(2m\pi/2)}{2m\pi}\,a_{2m}-2\bigl(a_{2m}^{2}+b_{2m}^{2}\bigr).
$$
Since $\sin(m\pi)=0$ for $m\in\mathbb Z$, the first term vanishes, so
$A_{2m}=-2(a_{2m}^{2}+b_{2m}^{2})$. By formulas (3.3) and (3.4) of Lemma 3
(the even-index branch), $a_{2m}=\tfrac12 c_{m}$ and $b_{2m}=\tfrac12 d_{m}$.
Therefore
$$
A_{2m}\;=\;-2\!\left(\tfrac{c_{m}^{2}}{4}+\tfrac{d_{m}^{2}}{4}\right)\;=\;-\tfrac12(c_{m}^{2}+d_{m}^{2}),
\qquad
A_{2m}^{2}\;=\;\tfrac14(c_{m}^{2}+d_{m}^{2})^{2}. \tag{5}
$$
Substituting (5) into (4):
$$
\frac14\sum_{m=1}^{T}(c_{m}^{2}+d_{m}^{2})^{2}\;\le\;\frac{\Omega-1/4}{2}
\;\Longleftrightarrow\;
\sum_{m=1}^{T}(c_{m}^{2}+d_{m}^{2})^{2}\;\le\;2\Omega-\tfrac12.
$$
This is exactly $(\star)$. $\blacksquare$

## 4. Convexity of $(\star)$

For each $m$ the function $(c_{m},d_{m})\mapsto (c_{m}^{2}+d_{m}^{2})^{2}$ is
the squared Euclidean norm raised to the fourth power, equivalently
$\lVert(c_{m},d_{m})\rVert_{2}^{4}$. The map $t\mapsto t^{4}$ on $t\ge0$ is
convex and non‑decreasing, and the Euclidean norm is convex; composition of
a non‑decreasing convex function with a convex function is convex. Hence
$(c_{m}^{2}+d_{m}^{2})^{2}$ is convex in $(c_{m},d_{m})$. Summing convex
functions preserves convexity, and the right-hand side $2\Omega-\tfrac12$
is affine. Therefore $(\star)$ defines a **convex** feasible region in the
variables $(c_{1},d_{1},\dots,c_{T},d_{T},\Omega)$.

## 5. Why this is a genuinely new constraint

The existing program does encode $\sum (c_{k}^{2}+d_{k}^{2})\le 1/2$
(Parseval on $[-1,1]$ for $f$, constraint (5.11)). Constraint $(\star)$ is
**not** implied by it: $(\star)$ involves the $\ell^{4}$‑norm of the vector
$\bigl((c_{m}^{2}+d_{m}^{2})\bigr)_{m\ge1}$ and *also* ties it to $\Omega$.
Concretely:

* The known constraint only bounds the $\ell^{1}$-norm of
  $(c_{m}^{2}+d_{m}^{2})_{m}$, allowing arbitrary concentration into a
  single coordinate.
* $(\star)$ bounds the $\ell^{2}$-norm of the same sequence — penalising
  concentration — and the bound shrinks linearly with $\Omega$.

Equivalently, $(\star)$ says
$\sum A_{2m}^{2}\le (\Omega-1/4)/2$: every individual even cosine
coefficient of $M$ — and the whole sequence in $\ell^{2}$ — is forced to
be small whenever $\Omega$ is close to $1/4$. This is information of a
qualitatively new type compared with the (linear) sign constraint
$A_{2m}\le 0$ used in (5.5).

## 6. Implementation outline

We introduce $T$ auxiliary non-negative variables
$s_{m}\ge c_{m}^{2}+d_{m}^{2}$ (each a single rotated SOC constraint that
CVXPY recognises as DCP because the right-hand side is a sum of squares
of affine variables). Then $(\star)$ becomes
$$
\sum_{m=1}^{T} s_{m}^{2}\;\le\;2\Omega-\tfrac12,
$$
which is a single second‑order‑cone constraint in $(s_{1},\dots,s_{T},\Omega)$.

In short: the addition is *one* quadratic-over-affine
constraint plus $T$ scalar slack epigraph constraints, both compatible
with MOSEK / SCS and matching the SOCP formulation White uses in
Appendix I.
