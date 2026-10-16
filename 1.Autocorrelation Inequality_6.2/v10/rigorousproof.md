# v10 — Rigorous proof and observed failure

## Setup

v10 attempts to extend v9's Fejér-kernel test to **non-even** $f$ by
introducing both $a_k$ (cosine) and $b_k$ (sine) Fourier coefficients
of $f$, plus a moment matrix $M$ on $(1, a_1, \dots, a_K, b_1, \dots, b_K)$
of size $2K + 1$.

## The Fejér constraint at multiple test points

For any test point $t_* \in [-1/2, 1/2]$ and the Fejér kernel $F_K$
(which is nonneg with $\int F_K = 1$):

$$
(F_K * F)(t_*) \;\le\; \|F_K\|_1 \cdot \|F\|_\infty \;=\; \Omega.
$$

Computing $(F_K * F)(t_*)$ via Plancherel:
$$
(F_K * F)(t_*) \;=\; 1 + 2\sum_{k=1}^{K}\Bigl(1 - \frac{k}{K+1}\Bigr)
\Bigl[(a_k^2 - b_k^2)\cos(2\pi k t_*) + 2 a_k b_k \sin(2\pi k t_*)\Bigr].
$$

Lifting via the moment matrix $M$ with $M[k,k] = a_k^2$,
$M[K+k, K+k] = b_k^2$, $M[k, K+k] = a_k b_k$ (in the rank-one truth):

$$
\Omega \;\ge\; 1 + 2\sum_{k=1}^{K}\Bigl(1 - \frac{k}{K+1}\Bigr)
\Bigl[(M[k,k] - M[K+k, K+k])\cos(2\pi k t_*)
+ 2 M[k, K+k] \sin(2\pi k t_*)\Bigr].
$$

For each test point $t_*$, this is a linear constraint on $(M, \Omega)$.
The rank-one $M^* = vv^{\top}$ from any admissible $f$ satisfies it
with the value being exactly $(F_K * F)(t_*) \le \Omega^* = \|F\|_\infty$.

So the constraints are **valid** for the relaxation.

## Why v10 fails numerically

At $N = 8, 16$ and $K = 8, 16, n_\text{tests} = 33$, MOSEK returns
$\Omega^{(10)} = 1$.

The reason is that the SDP lift $M$ does **not** constrain the diagonal
$M[K+k, K+k]$ (the lifted $b_k^2$) from above.  In particular, the
PSD lift only gives $M[K+k, K+k] \ge b_k^2$, so the relaxation can
freely **inflate** $M[K+k, K+k]$ to make the difference
$(M[k,k] - M[K+k, K+k])$ arbitrarily negative whenever
$\cos(2\pi k t_*) > 0$, and arbitrarily positive whenever
$\cos(2\pi k t_*) < 0$.  Combined with the cross term
$2 M[k, K+k]$ which can be inflated similarly via the off-diagonal
PSD freedom, the relaxation can make the right-hand side be exactly
$1$ at every test point, giving $\Omega = 1$.

This is the **same Fourier-side inflation degeneracy** noted in the
v2 verdict.  The v9 fix was to assume $f$ even (which forces
$b_k = 0$ and pins $M[K+k, K+k] = 0$).  v10 needs a similar
"upper bound on $M[K+k, K+k]$" mechanism for the non-even case,
which appears to require adding constraints from a different
direction.

## Validity and conclusion

The proposed constraints in v10 are **mathematically valid** (they
hold for the rank-one lift of any admissible $f$).  However, the
relaxation as written does not exclude the inflated-$M$ degeneracy,
so the LP optimum is the trivial $\Omega = 1$.

v10 is therefore "rigorous but vacuous": it neither violates the
relaxation principle nor improves the bound.  v11 will replace the
SDP lift with a tighter alternative — specifically, the Bochner /
Toeplitz PSD constraint on $(a_k, b_k)$, which ensures that the
sequence $(\hat f(k))_k$ is the Fourier sequence of an actual
non-negative function and thereby caps the off-diagonal moment
inflation.
