# v14 — Rigorous proof: autocorrelation Fejér constraint, valid for ALL admissible f

## Goal

v14 produces the **first rigorous lower bound on $C_{6.2}$ in this
project that does not require an even-$f$ restriction**.  The
construction follows the v9 verdict §5 recommendation: use the
Fejér-kernel test functional applied to the **autocorrelation**
$F * \check F$ where $\check F(t) = F(-t)$, instead of $F$ itself.
The autocorrelation has manifestly nonneg Fourier-side data, which
side-steps the M[K+k, K+k] inflation degeneracy that defeated v10–v13.

## Setup

* WLOG (justified in v7 §2): $f \ge 0$,
  $\operatorname{supp} f \subset [-1/4, 1/4]$, $\int f = 1$.
* No even restriction.
* $F = f * f$, $\Omega = \|F\|_\infty$.

## The autocorrelation Fejér test

**Definition.**  $\check F(t) := F(-t)$, the time-reversed $F$.
$F * \check F$ is the **autocorrelation** of $F$.

**Three properties.**

(A) $(F * \check F)\widehat{}(k) = \hat F(k) \hat F(-k)
   = \hat F(k) \overline{\hat F(k)} = |\hat F(k)|^{2}$.  This is
   real and nonneg for every $k$.

(B) $F * \check F \ge 0$ pointwise (a convolution of two nonneg
    functions, since both $F$ and $\check F$ are nonneg).

(C) $\|F * \check F\|_\infty = (F * \check F)(0) = \int F^{2}\,dt$.
   (The autocorrelation peaks at $0$.)

**The constraint chain.**

$$
\int F_K(t)\,(F * \check F)(t)\,dt
\;\overset{\text{(C')}}{\le}\;
\|F_K\|_1 \cdot \|F * \check F\|_\infty
\;\overset{\text{(C)}}{=}\;
\int F^{2}\,dt
\;\overset{\text{(D)}}{\le}\;
\Omega \cdot \int F\,dt
\;=\; \Omega.
$$

* (C') is Hölder $\int gh \le \|g\|_1\|h\|_\infty$.  $\|F_K\|_1 = 1$
  by (F2) of v9.
* (D) is $\int F^{2} = \int F\cdot F \le \|F\|_\infty \cdot \int F$
  using $F \ge 0$, $F \le \Omega$, $\int F = 1$.

The left-hand side, by Parseval, is
$$
\int F_K (F * \check F)\,dt
\;=\; \sum_{|k|\le K}\widehat{F_K}(k)\,\widehat{(F*\check F)}(k)
\;=\; \sum_{|k|\le K}\Bigl(1 - \frac{|k|}{K+1}\Bigr)\,|\hat F(k)|^{2}.
$$
Since $|\hat F(k)|^{2} = |\hat F(-k)|^{2}$ and $|\hat F(0)|^{2} = 1$,
the sum simplifies to
$$
1 + 2\sum_{k=1}^{K}\Bigl(1 - \frac{k}{K+1}\Bigr)\,|\hat F(k)|^{2}.
$$

Combining,
$$
1 + 2\sum_{k=1}^{K}\Bigl(1 - \frac{k}{K+1}\Bigr)\,|\hat F(k)|^{2}
\;\le\; \Omega.
\tag{$\star$}
$$

**Key feature.**  The Fejér weights are nonnegative, and so are the
$|\hat F(k)|^{2}$ — every term in $(\star)$ has the **same sign**.
There is no inflation degeneracy.

## Lifting to a convex constraint

We must express $|\hat F(k)|^{2}$ in terms of LP variables.  Recall
$\hat F(k) = \hat f(k)^{2} = (a_k - i b_k)^{2}$, so
$$
|\hat F(k)|^{2}
\;=\; |\hat f(k)^{2}|^{2}
\;=\; |\hat f(k)|^{4}
\;=\; (a_k^{2} + b_k^{2})^{2}.
$$

Introduce two layers of LP variables.

**Layer 1.**  Let $u_k := M[k,k] + M[K+k, K+k]$ (linear in the SDP-lift
$M$).  In the rank-one truth, $M[k,k] = a_k^{2}$ and
$M[K+k, K+k] = b_k^{2}$, so $u_k = a_k^{2} + b_k^{2} = |\hat f(k)|^{2}$.

In the relaxation, the PSD lift on $M$ gives only
$M[k,k] \ge a_k^{2}$ and $M[K+k, K+k] \ge b_k^{2}$, so
$u_k \ge a_k^{2} + b_k^{2}$.  This is the *correct* direction for
us: the LP wants to **minimize** $u_k$ to make $\Omega$ small, so it
will set $u_k = a_k^{2} + b_k^{2}$ at the rank-one boundary.

**Layer 2.**  Let $v_k \ge u_k^{2}$, a rotated SOC (convex)
inequality.  In the rank-one truth, $v_k = u_k^{2} = (a_k^{2} + b_k^{2})^{2} = |\hat F(k)|^{2}$.

Since the LP minimizes $\Omega$, and $v_k$ appears with positive
weight in the constraint $(\star\star)$ below, the LP sets $v_k$ to
its smallest value satisfying the SOC, i.e. $v_k = u_k^{2}$.

**Lifted constraint.**

$$
\Omega \;\ge\; 1 + 2\sum_{k=1}^{K}\Bigl(1 - \frac{k}{K+1}\Bigr)\,v_k.
\tag{$\star\star$}
$$

## Validity

Every admissible $f$ — *not just even* — produces a feasible
$(p, a, b, M, u, v, \Omega^{*} = \|f * f\|_\infty)$:

* Rank-one $M = vv^{\!\top}$ with $v = (1, a_1, \dots, a_K, b_1, \dots, b_K)$;
* $u_k = a_k^{2} + b_k^{2}$ at the PSD-lift boundary;
* $v_k = u_k^{2} = |\hat F(k)|^{2}$ at the SOC boundary;
* $\Omega^{*}$ satisfies $(\star\star)$ at equality, by the
  derivation of $(\star)$ above.

So the v14 optimum $\Omega^{(14)} \le \Omega^{*}$ for every admissible
$f$, giving
$$
\Omega^{(14)} \;\le\; C_{6.2}.
$$

This is a **rigorous lower bound on $C_{6.2}$** (the unrestricted
constant), the first one in the project.

## Numerical results

| $N$ | $K$ | $\Omega^{(14)}$ |
|-----|-----|-----------------|
| 8   | 4   | 1.026           |
| 8   | 8   | 1.033           |
| 8   | 16  | 1.036           |
| 16  | 4   | 1.040           |
| 16  | 8   | 1.055           |
| 16  | 16  | 1.063           |

Strictly above $1$ for every $N \ge 8$.  The bound grows with both
$N$ and $K$, and should converge to a non-trivial limit (less than
2) as both grow.

The bound is well below the SOTA $C_{6.2} \ge 1.28$ at the small
$(N, K)$ tested here.  Larger $N, K$ runs are needed to see how close
v14 can get; the next iteration (v15) will explore larger
parameters and add the v11 windowed (W) constraints to tighten.

## What's special about the autocorrelation lift

The conventional Fejér test on $F$ alone gives a constraint of the
form $\Omega \ge 1 + 2\sum (\text{weight}) \cdot (a_k^{2} - b_k^{2})$,
where the *signed* difference $a_k^{2} - b_k^{2}$ admits inflation
of the $b_k^{2}$ term (via $M[K+k, K+k]$ in the lift).  This
defeated v2, v8, v10, v12, v13.

The autocorrelation Fejér test gives a constraint of the form
$\Omega \ge 1 + 2\sum (\text{weight}) \cdot |\hat F(k)|^{2}$, where
$|\hat F(k)|^{2} = (a_k^{2} + b_k^{2})^{2}$ is *intrinsically nonneg*.
The "minus" sign that allowed the inflation degeneracy is gone.

This is the **right structural fix** that the v9 verdict pointed to.
