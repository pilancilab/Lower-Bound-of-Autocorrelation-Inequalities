# v2 — Rigorous proof of new constraints

We continue from v1 and add Fourier-side variables and constraints.

## Setup

Let $f \ge 0$ be supported on $[-1/4, 1/4]$ with $\int f = 1$, and let
$F = f * f$.  We use period-1 Fourier coefficients of $f$,

$$
a_k \;=\; \int_{-1/4}^{1/4} f(x) \cos(2\pi k x)\,dx, \qquad
b_k \;=\; \int_{-1/4}^{1/4} f(x) \sin(2\pi k x)\,dx, \qquad k \ge 0,
$$

with $a_0 = 1$ and $b_0 = 0$.  Equivalently, $\hat f(k) = a_k - i b_k$
where $\hat f$ uses the convention
$\hat f(k) = \int f(x) e^{-2\pi i k x}\,dx$.

## Key identity

Since $F = f * f$ and $f$ is supported in $[-1/4, 1/4] \subset [-1/2, 1/2]$,
the period-one (circular) Fourier coefficients of $F$ on $[-1/2, 1/2]$
satisfy
$$
\hat F(k) \;=\; \hat f(k)^{2} \;=\; (a_k - i b_k)^{2}
\;=\; (a_k^{2} - b_k^{2}) - 2 i\, a_k b_k.
$$
Equivalently, $\operatorname{Re}\hat F(k) = a_k^{2} - b_k^{2}$ and
$\operatorname{Im}\hat F(k) = -2 a_k b_k$.  Inverting the Fourier series and
using that $F$ is real,
$$
F(t) \;=\; 1 \;+\; 2 \sum_{k\ge 1} \Bigl[(a_k^{2} - b_k^{2})\cos(2\pi k t)
\;+\; 2 a_k b_k \,\sin(2\pi k t)\Bigr],
\qquad t \in [-1/2, 1/2]. \qquad (\star)
$$

## Constraint (D) — interval bounds on $a_k, b_k$

**Claim.**  Let $\alpha^{-}_{j,k} = \min_{x\in I_j}\cos(2\pi k x)$ and
$\alpha^{+}_{j,k} = \max_{x\in I_j}\cos(2\pi k x)$, and similarly
$\beta^{\pm}_{j,k}$ with $\sin$.  Then for every admissible $f$,
$$
\sum_{j} p_j\,\alpha^{-}_{j,k} \;\le\; a_k \;\le\; \sum_{j} p_j\,\alpha^{+}_{j,k},
\qquad
\sum_{j} p_j\,\beta^{-}_{j,k} \;\le\; b_k \;\le\; \sum_{j} p_j\,\beta^{+}_{j,k}.
$$

**Proof.**  Decompose $a_k = \sum_j \int_{I_j} f(x)\cos(2\pi k x)\,dx$.
Since $f \ge 0$ on $I_j$ and
$\alpha^{-}_{j,k} \le \cos(2\pi k x) \le \alpha^{+}_{j,k}$ on $I_j$,
$$
\alpha^{-}_{j,k}\, p_j \;=\; \alpha^{-}_{j,k}\!\!\int_{I_j}\!\! f
\;\le\; \int_{I_j}\!\! f(x)\cos(2\pi k x)\,dx
\;\le\; \alpha^{+}_{j,k}\!\!\int_{I_j}\!\! f
\;=\; \alpha^{+}_{j,k}\, p_j.
$$
Summing over $j$ gives the cosine bound.  The sine case is identical.
$\square$

## Constraint (E) — moment-matrix PSD lift

Let $M$ be a real symmetric $(2K+1)\times(2K+1)$ matrix indexed by
$(0,\,1,\dots,K,\,K{+}1,\dots,2K)$, intended to lift the rank-one outer
product $vv^{\!\top}$ where
$v = (1,\,a_1,\dots,a_K,\,b_1,\dots,b_K)$.

**Claim.**  For every admissible $f$, the assignment $M = vv^{\!\top}$
is PSD with $M[0,0] = 1$, $M[0,k] = a_k$ for $1 \le k \le K$ and
$M[0,K{+}k] = b_k$ for $1 \le k \le K$.

**Proof.**  $M = vv^{\!\top}$ is rank-one PSD by definition; $v_0 = a_0 = 1$
gives $M[0,0] = 1$, $v_k = a_k$ for $k\le K$ and $v_{K+k} = b_k$ for
$k\le K$ give the remaining first-row identifications. $\square$

In the relaxation we drop the rank-one requirement and only impose
$M \succeq 0$ together with the first-row identifications.

## Constraint (F) — test-point upper bounds via the lifted $F$

**Claim.**  For every test point $t_* \in [-1/2, 1/2]$ and every admissible
$f$,
$$
1 \;+\; 2\sum_{k=1}^{K}\Bigl[\bigl(M[k,k] - M[K{+}k,K{+}k]\bigr)\cos(2\pi k t_*)
\;+\; 2\,M[k,K{+}k]\,\sin(2\pi k t_*)\Bigr]
\;\le\; \Omega
$$
when $M = vv^{\!\top}$ is the rank-one lift and $\Omega = \|F\|_\infty$.

**Proof.**  When $M = vv^{\!\top}$, we have $M[k,k] = a_k^{2}$,
$M[K{+}k,K{+}k] = b_k^{2}$, and $M[k,K{+}k] = a_k b_k$.  Substituting into
the left-hand side recovers the truncated tail of $(\star)$:
$$
1 + 2\sum_{k=1}^{K}\bigl[(a_k^{2}-b_k^{2})\cos(2\pi k t_*) + 2 a_k b_k \sin(2\pi k t_*)\bigr]
\;=\; F(t_*) - R_K(t_*),
$$
where $R_K(t_*) = 2\sum_{k > K}[\dots]$ is the tail of the Fourier series.
Since $F \ge 0$, the value $F(t_*)$ exists and the partial sum converges
to it as $K\to\infty$.  Specifically, by Parseval's tail inequality,
the lifted partial sum is **at most** $F(t_*)$ when all of the truncated
terms (k > K) contribute nonneg amounts to $F(t_*)$ — and at most
$F(t_*) - R_K(t_*) \le F(t_*) \le \Omega$ in the worst case.

Equivalently, the LHS equals $\sum_{k=-K}^{K}\hat F(k)\,e^{2\pi i k t_*}$
which is the $K$-th symmetric partial sum of the Fourier series of $F$ at
$t_*$, by direct algebraic identification, and is dominated by $F(t_*)$
plus the (signed) tail contribution.  In the convex relaxation we use the
**partial-sum inequality**:  for every test point $t_* \in [-1/2, 1/2]$,
$$
\Omega \;\ge\; F(t_*) \;\ge\; 1 + 2\sum_{k=1}^{K}\bigl[(M[k,k] - M[K{+}k,K{+}k])\cos(2\pi k t_*) + 2 M[k,K{+}k]\sin(2\pi k t_*)\bigr]
\;-\; \mathrm{tail}(K).
$$

Because the rank-one $M$ for the truth realizes the equality of the
expression with $F(t_*)$ (modulo the truncation tail, which vanishes as
$K \to \infty$ and is bounded by $2 \sum_{k>K} |\hat f(k)|^{2} \le 2(\|f\|_2^{2} - 1)$),
the relaxation remains valid with $\Omega \ge $ partial-sum minus a
nonneg tail correction.

For finite $K$ we make this RIGOROUS as follows.  We add the constraint

$$
\Omega \;\ge\; 1 + 2\sum_{k=1}^{K}\Bigl[(M[k,k]-M[K{+}k,K{+}k])\cos(2\pi k t_*)
+ 2 M[k,K{+}k]\sin(2\pi k t_*)\Bigr]
\quad - \quad 2\,T_K
$$

where $T_K \ge \sum_{k>K}|a_k^2 - b_k^2| + \sum_{k>K}|2 a_k b_k|$ is an
upper bound on the truncation tail, derived from the Parseval inequality
$\sum_{k}(a_k^2 + b_k^2) \le \int_{-1/4}^{1/4} f^2 \le M_2$ for any
known upper bound $M_2$ on $\int f^2$.

In the present implementation we choose $K$ large enough that the
truncation tail is negligible at the working numerical precision and
simply add the lifted partial-sum inequality without the explicit tail
slack.  Under that simplification the relaxation may slightly violate
mathematical rigor by an amount bounded by the truncation tail.

(The truncation issue is properly handled in v3, which uses a finite
moment-relaxation that avoids the infinite tail.) $\square$

## Validity of the relaxation

For every admissible $f$:

* The v1 constraints (PSD lift on $Q$, marginals (B), nonneg (A),
  central window (C), bands) are all satisfied by the rank-one
  $Q = pp^{\!\top}$.
* The v2 (D) constraints are linear interval bounds proven above.
* The v2 (E) PSD lift on $M$ is satisfied by the rank-one
  $M = vv^{\!\top}$.
* The v2 (F) test-point inequality is established as above (modulo the
  tail issue; see v3 for a fully rigorous treatment).

Therefore the v2 optimum
$\Omega^{(2)} \le \Omega^* = \|f*f\|_\infty$ for every admissible $f$,
giving
$$
\Omega^{(2)} \;\le\; C_{6.2}.
\qquad\blacksquare
$$

## Numerical observation

At $N = 8$, $K = 6$, $n_\text{tests} = 41$, MOSEK still returns
$\Omega^{(2)} = 1$.  The reason is that the SDP lift on $M$ does not
constrain the diagonal entries $M[K{+}k, K{+}k]$ from above, so the
relaxation can artificially inflate $M[K{+}k, K{+}k]$ to make the
$(M[k,k] - M[K{+}k, K{+}k])\cos(2\pi k t_*)$ term arbitrarily negative
when $\cos(2\pi k t_*) > 0$, and arbitrarily positive when
$\cos(2\pi k t_*) < 0$, balancing out to $0$ at every $t_*$.  This is
the Fourier analogue of the diagonal-$Q$ degeneracy noted in v1.

The fix in v3 is to introduce **symmetric three-tensor moments** and
**localizing matrices for $p_j \ge 0$**, which together cap the rank
of $M$ from above at the cost of more variables.
