# v12 — Rigorous proof and observed failure

## Goal

v12 attempts the v7 verdict's recommendation #2: build the F-side
LP with $W_m$ as cell averages of $F$ and link to f-side via the
rank-one identification $\hat F(k) = \hat f(k)^{2}$, lifted via the
PSD lift on $(a_k, b_k)$.

## Constraints (rigorous)

(P0) Mass: $\sum p = 1$, $L_F \sum W = 1$, $W_m \le \Omega$.
(P1) PSD lift on $M$: $M \succeq 0$, $M[0,0] = 1$,
     $M[0, k] = a_k$, $M[0, K+k] = b_k$.
(P2) f-side bounds on $a_k$, $b_k$ from $p$ (linear, via min/max
     of cos/sin on each $I_j$).
(P3) F-side bounds on $A_k = a_k^2 - b_k^2$ (lifted as
     $M[k,k] - M[K+k, K+k]$) and $B_k = -2 a_k b_k$ (lifted as
     $-2 M[k, K+k]$) from $W$ (linear, via min/max of cos/sin on $B_m$):
     $L_F \sum_m W_m (\min\cos)_{m,k} \le M[k,k] - M[K+k, K+k] \le L_F \sum_m W_m (\max\cos)_{m,k}$
     and similarly for $B_k$.

All constraints (P0)–(P3) are rigorously valid: the rank-one lift
$M^* = vv^{\top}$ from any admissible $f$ satisfies them.

## Why v12 fails

At $N = 8, 12$ and $K = 6, 12$, MOSEK returns $\Omega^{(12)} = 1$
(the trivial bound).

The reason is the **M[K+k, K+k] inflation degeneracy** that defeated
v10:

* The PSD lift on $M$ only gives $M[K+k, K+k] \ge b_k^2$ (lower bound).
* The relaxation can freely set $M[K+k, K+k]$ to any value satisfying
  $M[K+k, K+k] \ge b_k^2$ AND the F-side bound from $W$.
* By inflating $M[K+k, K+k]$, the lifted $A_k = M[k,k] - M[K+k, K+k]$
  becomes arbitrarily small (or negative), so the F-side constraint
  $A_k \ge L_F \sum_m W_m (\min \cos)_{m,k}$ is satisfied for *uniform*
  $W = 1$ even though the truth would force $W$ to be non-uniform.

The result is that the relaxation finds $(p, W, M)$ with $\Omega = 1$
and $W$ uniform, all rigorously satisfying the constraints.

The fix requires **bounding $M[K+k, K+k]$ from above** by something
that is small for non-even $f$.  Two options:

1. **Even-restriction**: enforce $b_k = 0$ and $M[K+k, K+k] = 0$.
   This is what v9 does and gives a non-trivial bound on
   $C^{\mathrm{even}}_{6.2}$, but does not bound $C_{6.2}$.

2. **Bochner Toeplitz lift on (a_k, b_k)**: add the constraint that
   the Hermitian Toeplitz matrix $[\hat f(j-k)]_{j,k}$ is PSD (with
   $\hat f(k) = a_k - i b_k$, $\hat f(-k) = \hat f(k)^*$, $\hat f(0) = 1$).
   This is the Bochner condition encoding $f \ge 0$.  In the relaxation,
   it pins down $(a_k, b_k)$ to genuine Fourier sequences of nonneg
   measures, which limits the inflation freedom on $M$.

v13 will add the Bochner Toeplitz constraint.

## Validity

Despite failing numerically, the v12 constraints are **mathematically
valid**: every admissible $f$ produces a feasible $(p, W, M)$ at
$\Omega^* = \|f * f\|_\infty$.  So v12 is rigorous but vacuous.
