# v3 — Rigorous proof of the Lasserre level-2 lift

## Setup

We carry forward the v0 setting and notation: $f \ge 0$, $\operatorname{supp} f \subset [-1/4, 1/4]$,
$\int f = 1$, and we discretize $f$ by setting $p_j = \int_{I_j} f$ on the
$2N$ intervals of width $L = 1/(4N)$ tiling $[-1/4,1/4]$.  We minimize
$\Omega \ge \|f * f\|_\infty$.

The variables of v3 are real numbers $y_\alpha$ indexed by multi-indices
$\alpha \in \mathbb N^{2N}$ with $|\alpha| \le 4$, intended to lift the
moments

$$
y_\alpha \;:=\; \int p^\alpha\, d\mu(p) \;=\; \prod_j p_j^{\alpha_j}
\quad\text{(in the rank-one truth at $p$).}
$$

In particular, $y_{\mathbf 0} = 1$, $y_{e_j} = p_j$, and
$y_{e_j + e_k} = p_j p_k$ — these recover the variables of v0/v1/v2.

## Constraints (added relative to v0)

### (G) Moment-matrix PSD lift at level 2

Let $\alpha^{(1)}, \dots, \alpha^{(s_2)}$ enumerate the multi-indices of
degree $\le 2$.  For $n = 2N$ this gives $s_2 = 1 + n + \binom{n+1}{2}$
indices ($s_2 = 28$ at $n = 6$).  Define the moment matrix

$$
M_2[i, j] \;:=\; y_{\alpha^{(i)} + \alpha^{(j)}}, \qquad 1 \le i, j \le s_2.
$$

**Claim.**  In the rank-one truth, $M_2 = vv^{\!\top}$ where
$v_i = p^{\alpha^{(i)}}$.  In particular $M_2 \succeq 0$.

**Proof.**  Direct computation: $(vv^{\!\top})[i,j] = v_i v_j = p^{\alpha^{(i)}}\!\cdot p^{\alpha^{(j)}} = p^{\alpha^{(i)} + \alpha^{(j)}} = y_{\alpha^{(i)}+\alpha^{(j)}}$.
Rank-one outer products are PSD. $\square$

In the relaxation we therefore impose $M_2 \succeq 0$ and drop the rank-one
requirement.

### (H) Localizing matrices for $p_j \ge 0$

For each $j \in \{1, \dots, n\}$ and the multi-indices $\alpha^{(1)}, \dots, \alpha^{(s_1)}$
of degree $\le 1$ (so $s_1 = 1 + n$), define

$$
L^{(j)}[i, k] \;:=\; y_{\alpha^{(i)} + \alpha^{(k)} + e_j}.
$$

**Claim.**  In the rank-one truth, $L^{(j)} = p_j\, w_j w_j^{\!\top}$ where
$w_j$ is the rank-one moment vector (degree $\le 1$).  In particular, since
$p_j \ge 0$, $L^{(j)} \succeq 0$.

**Proof.**  $L^{(j)}[i, k] = p^{\alpha^{(i)} + \alpha^{(k)} + e_j} = p_j \cdot p^{\alpha^{(i)}+\alpha^{(k)}} = p_j \cdot w_{j,i} w_{j,k}$
with $w_{j,i} = p^{\alpha^{(i)}}$ (which does not depend on $j$ at the
level-1 monomials except through scaling).  More precisely the truth
factors as $L^{(j)} = p_j \cdot (w w^{\!\top})$ for the level-1 moment
vector $w$, and $p_j (w w^{\!\top}) \succeq 0$ since $p_j \ge 0$ and
$ww^{\!\top}$ is rank-one PSD. $\square$

### (I) Linear marginal/equality constraints

The mass constraint $\sum_j p_j = 1$ propagates to all levels: for each
$\alpha$ with $|\alpha| \le 3$,

$$
\sum_{j=1}^{n} y_{\alpha + e_j} \;=\; y_\alpha.
$$

**Claim.**  In the rank-one truth this is exact.

**Proof.**  $\sum_j y_{\alpha + e_j} = \sum_j p_j \cdot p^\alpha = p^\alpha \sum_j p_j = p^\alpha \cdot 1 = y_\alpha$. $\square$

We impose this as a linear equality on the moment variables.

## Validity of the relaxation

For every admissible $f$, set $p_j^* = \int_{I_j} f$ and
$y_\alpha^* = \prod_j (p_j^*)^{\alpha_j}$.  Then:

* $y_{\mathbf 0}^* = 1$ and $\sum_j y_{e_j}^* = 1$ (mass).
* $M_2 = vv^{\!\top}$ with $v_i = (p^*)^{\alpha^{(i)}}$, hence $M_2 \succeq 0$.
* Each $L^{(j)} = p_j^* (ww^{\!\top}) \succeq 0$.
* The marginal/equality constraints (I) hold by direct computation.
* The v0 band-integral constraints
  $\tfrac12 \sum_{j+k\in\{m,m+1\}} y_{e_j + e_k} \le L \cdot \Omega$
  are inherited from v0 (with $y_{e_j+e_k} = p_j^* p_k^*$).

So the (Omega, y) coming from any admissible $f$ is feasible for v3, and
the v3 optimum $\Omega^{(3)} \le \Omega^* = \|f * f\|_\infty$.  Hence
$\Omega^{(3)} \le C_{6.2}$.

$\blacksquare$

## Numerical observation

At $N = 3$ ($n = 6$, with $210$ moments and $M_2$ of size $28\times 28$),
MOSEK still returns $\Omega^{(3)} = 1$.

The reason is structural and important.  Consider the probability measure
$\mu$ on the simplex $\Delta_n$ that is uniform on the basis vectors
$\{e_1, e_2, \dots, e_n\}$.  Under this measure,

$$
\mathbb E_\mu[p^\alpha] \;=\; \frac{1}{n} \cdot \#\{l \,:\, \operatorname{supp}\alpha \subset \{l\}\}
\;=\; \begin{cases} 1 & \alpha = \mathbf 0,\\ 1/n & \alpha = k e_l, \;\; k \ge 1, \; l \in [n],\\ 0 & |\operatorname{supp}\alpha| \ge 2.\end{cases}
$$

This moment sequence is feasible for v3 because:

* It is the moment sequence of a *bona fide* probability measure on
  $\Delta_n$, so all PSD constraints (the moment matrix and every
  localizing matrix) are automatically satisfied.
* The mass and marginal constraints hold by construction.
* The induced $Q[j,k] = \mathbb E[p_j p_k] = (1/n)\delta_{jk}$ is a
  diagonal matrix, and the band integrals reduce to $1/(2N) + 0 = 1/(2N)$
  on alternating bands, exactly matching $L\,\Omega$ with $\Omega = 1$.

This is the **same** degeneracy that defeated v1 and v2: the convex
relaxation cannot prevent the optimization from "averaging" the problem
over a distribution of $f$'s with very different $\|f * f\|_\infty$.

### Why higher Lasserre levels do *not* fix the issue

The basis-vector measure $\mu$ is a *real* probability measure on
$\Delta_n$.  Therefore its moment sequence $y_\alpha = \mathbb E_\mu[p^\alpha]$
is a valid moment sequence at every level $\ell$ — the level-$\ell$
moment matrix is the rank-$n$ Gram matrix of the basis-vector evaluations,
which is always PSD.  No finite Lasserre level can rule this out.

The key fix must therefore decouple $\Omega$ from a single $f$, e.g. by
introducing constraints that involve $\Omega$ together with the moments
in a way that no convex combination of point-mass distributions can
satisfy.  v4 addresses this by working entirely on the $F$-side via a
discrete LP that bounds the autocorrelation values directly using exact
band integrals computed from a Fourier-cosine representation, mirroring
White's Section 4 LP.
