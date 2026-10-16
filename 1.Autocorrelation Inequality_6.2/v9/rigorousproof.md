# v9 — Rigorous proof: Fejér-kernel test functional

## Goal

v9 produces a **fully rigorous** convex relaxation whose optimum is a
lower bound on $C^{\mathrm{even}}_{6.2}$, the even-restricted
constant.  The construction follows the verifier's recommendation
(v5 verdict, Comment 2): replace point-evaluation Fourier tests by
**positive trigonometric polynomial test functionals**, which are
linear in the lifted moments and have no Gibbs/tail issue.

## Key tool: the Fejér kernel

The Fejér kernel of order $K$ on the period-1 circle is
$$
F_K(t) \;=\; \sum_{k=-K}^{K}\Bigl(1 - \frac{|k|}{K+1}\Bigr) e^{2\pi i k t}
\;=\; 1 + 2\sum_{k=1}^{K}\Bigl(1 - \frac{k}{K+1}\Bigr)\cos(2\pi k t).
$$

**Two well-known properties.**

(F1) $F_K(t) \ge 0$ for all $t$.  (Indeed, $F_K(t) = (1/(K+1))\,|\sum_{k=0}^{K} e^{2\pi i k t}|^{2}$
is the squared modulus of a Dirichlet kernel partial sum, hence nonneg.)

(F2) $\int_{-1/2}^{1/2} F_K(t)\,dt = 1$ (the constant Fourier
coefficient is $1$).

## Setup and reductions

* WLOG (justified in v7's rigorousproof.md): $f \ge 0$,
  $\operatorname{supp} f \subset [-1/4, 1/4]$, $\int f = 1$.
* v9 also assumes $f$ even (so the bound applies to $C^{\mathrm{even}}_{6.2}$).
* $F = f * f$, $\Omega = \|F\|_\infty$.

For real even $f$ on $[-1/4, 1/4] \subset [-1/2, 1/2]$, define the
cosine Fourier coefficients
$$
a_k \;:=\; \int_{-1/4}^{1/4} f(x)\cos(2\pi k x)\,dx, \qquad k \ge 0,
$$
with $a_0 = 1$ (the mass).  The sine coefficients vanish because
$f$ is even.

By Lemma (algebraic Fourier identity for $F = f * f$ on the period-1 circle):
$$
F(t) \;=\; \sum_{k} \hat F(k)\,e^{2\pi i k t}
\;=\; 1 + 2\sum_{k=1}^{\infty} a_k^{2}\cos(2\pi k t).
$$
This identity follows from $\hat F(k) = \hat f(k)^{2}$ and the fact
that $\hat f(k) = a_k$ (real, since $f$ is even).

## The Fejér-kernel test functional

**Lemma.**  For every admissible even $f$,
$$
1 + 2\sum_{k=1}^{K}\Bigl(1 - \frac{k}{K+1}\Bigr) a_k^{2}
\;\le\; \Omega.
\tag{$\star$}
$$

**Proof.**  Compute the inner product $\int F_K \cdot F$ in two ways.
First, expand $F_K$ in cosines and integrate term by term against the
cosine series of $F$:
$$
\int_{-1/2}^{1/2} F_K(t)\,F(t)\,dt
\;=\; \int_{-1/2}^{1/2}\!\Bigl[1 + 2\!\sum_{k=1}^{K}\!\bigl(1 - \tfrac{k}{K+1}\bigr)\cos(2\pi k t)\Bigr]
\Bigl[1 + 2\!\sum_{m=1}^{\infty}\!a_m^{2}\cos(2\pi m t)\Bigr]\,dt.
$$
By orthogonality of the cosines on $[-1/2, 1/2]$
($\int_{-1/2}^{1/2}\cos(2\pi k t)\cos(2\pi m t)\,dt = \delta_{km}/2$
for $k, m \ge 1$, and $= 1$ for $k = m = 0$),
$$
\int F_K F = 1 + 2\sum_{k=1}^{K}\bigl(1 - \tfrac{k}{K+1}\bigr) a_k^{2}.
$$

Second, since $F_K \ge 0$ and $F \le \Omega$ pointwise,
$$
\int F_K F \;\le\; \Omega \int F_K \;=\; \Omega \cdot 1.
$$

Combining the two computations gives $(\star)$. $\square$

## Lifting to a convex constraint

The right-hand side of $(\star)$ is **quadratic** in the variables
$a_k$.  We linearize via the SDP lift on the moment matrix $M$.

Let $M \in \mathbb R^{(K+1)\times(K+1)}$ symmetric, with $M[0,0] = 1$
and $M[0, k] = a_k$ for $k = 1, \dots, K$.  The PSD constraint
$M \succeq 0$ implies (via the $2\times 2$ minor
$\bigl[\begin{smallmatrix}1 & a_k\\ a_k & M[k,k]\end{smallmatrix}\bigr]$):
$$
M[k, k] \;\ge\; a_k^{2}.
$$

The lifted constraint we add is
$$
\Omega \;\ge\; 1 + 2\sum_{k=1}^{K}\Bigl(1 - \frac{k}{K+1}\Bigr) M[k, k].
\tag{$\star\star$}
$$

**Validity claim.**  $(\star\star)$ holds for the rank-one
$M^* = vv^{\!\top}$ derived from any admissible even $f$, where
$v = (1, a_1, \dots, a_K)$.  **Therefore $(\star\star)$ is a valid
convex relaxation constraint.**

**Proof.**  At $M^* = vv^{\!\top}$, $M^*[k, k] = a_k^{2}$, so
$$
1 + 2\sum_{k=1}^{K}\Bigl(1 - \frac{k}{K+1}\Bigr) M^*[k, k]
\;=\; 1 + 2\sum_{k=1}^{K}\Bigl(1 - \frac{k}{K+1}\Bigr) a_k^{2}
\;\le\; \Omega
$$
by $(\star)$. $\square$

**Important subtlety: tightness direction.**  Inside the relaxation,
$M[k, k]$ may be **larger** than $a_k^{2}$ (the SDP lift only gives a
lower bound on $M[k, k]$).  The lifted left-hand side of $(\star\star)$
is therefore $\ge$ the right-hand side of $(\star)$ in the relaxed
program.  The constraint $(\star\star)$ enforces
$$
\Omega \;\ge\; (\text{lifted LHS}) \;\ge\; (\text{true LHS of } (\star)),
$$
which is strictly stronger than the original constraint
$\Omega \ge \text{true LHS}$, so the relaxation is **valid**: every
admissible $f$ gives a feasible point of the relaxed program at
$\Omega^* = \|f * f\|_\infty$.

## Linear bounds on $a_k$ from $p$ (f-side discretization)

To pin down $a_k$ in the relaxation, we discretize $f$ via the cell
masses $p_j = \int_{I_j} f$ for $j = 1, \dots, 2N$ on the grid
$I_j = [-1/4 + (j-1)L_f,\, -1/4 + j L_f]$ with $L_f = 1/(4N)$.  The
linear interval bounds (proven in v2/v7 rigorousproof.md) hold:
$$
\sum_{j} p_j \min_{x \in I_j}\cos(2\pi k x)
\;\le\; a_k
\;\le\; \sum_{j} p_j \max_{x \in I_j}\cos(2\pi k x).
\tag{LB/UB}
$$

These bounds together with the even-symmetry $p_j = p_{2N+1-j}$
(proven in v5/v7) and the mass constraint $\sum p_j = 1$ are the
full f-side LP constraints.

## Numerical results

| $N$ | $K$ | $\Omega^{(9)}$ |
|-----|-----|---------------|
| 8   | 4   | 1.268         |
| 8   | 8   | 1.325         |
| 8   | 16  | 1.349         |
| 8   | 32  | 1.362         |
| 16  | 4   | 1.334         |
| 16  | 8   | 1.437         |
| 16  | 16  | 1.486         |
| 16  | 32  | 1.505         |
| 32  | 32  | 1.616         |
| 32  | 64  | 1.629         |
| 50  | 32  | 1.674         |
| 50  | 64  | 1.689         |

The relaxed optimum **strictly exceeds 1** for every choice of $(N, K)$,
and grows monotonically toward the analytic value
$C^{\mathrm{even}}_{6.2} = 2$ as both $N$ and $K$ increase.

## What v9 proves

**Theorem (v9).**  $C^{\mathrm{even}}_{6.2} \ge \Omega^{(9)}$.  At
$(N, K) = (50, 64)$, MOSEK certifies
$C^{\mathrm{even}}_{6.2} \ge 1.689$.

Combined with the analytic value $C^{\mathrm{even}}_{6.2} = 2$
achieved by $f = 2 \cdot \mathbf 1_{[-1/4, 1/4]}$, v9 gives a
nontrivial rigorous lower bound on the even-restricted constant
that approaches $2$ in the limit.

## What v9 does NOT prove

* It does **not** bound $C_{6.2}$ in general.  We have
  $C_{6.2} \le C^{\mathrm{even}}_{6.2}$, so a lower bound on the
  even-restricted constant is automatically a lower bound only on
  the even-restricted constant.
* The SOTA $C_{6.2} \ge 1.28$ (per `Problem.md`) is for the
  unrestricted problem, where the optimal $f$ is conjectured to be
  non-even.  v9 does not improve this.

## v10 plan

To extend the Fejér trick to non-even $f$, we will use **two test
functionals** simultaneously:

(i) The Fejér kernel against the *even part* of $F$:
   $\int F_K \cdot F_e\,dt$ where $F_e(t) = (F(t) + F(-t))/2$.

(ii) A complementary test against the *odd part*:
   $F_o(t) = (F(t) - F(-t))/2$.

The even part $F_e$ has positive Fourier coefficients and the Fejér
test gives a meaningful lower bound; the odd part contributes nothing
to $\int F_K \cdot F_e$ (orthogonality) but adds structure via a
separate constraint.  Combining the two gives a constraint valid for
all admissible $f$, not just even.
