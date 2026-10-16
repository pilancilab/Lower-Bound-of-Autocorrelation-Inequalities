# v47 — Trig Cauchy–Schwarz constraints relating $M$ to $a_{2k}$

## Goal

In v45 the moment matrix $M$ is constrained only by $M \succeq 0$
and by first-row equalities $M[0,k] = a_k$, $M[0,K+k] = b_k$.
This is the standard "level-1 Lasserre" relaxation.  The diagonal
entries $M[k,k]$ and $M[K+k,K+k]$ are SDP relaxations of $a_k^2$
and $b_k^2$, and the $M \succeq 0$ constraint forces $M[k,k] \ge a_k^2$
(by Schur on the $2 \times 2$ submatrix
$\bigl[\begin{smallmatrix} 1 & a_k \\ a_k & M[k,k] \end{smallmatrix}\bigr]$).

But this is the only constraint on $M[k,k]$ from above — and the
SDP is free to push $M[k,k]$ as high as needed to weaken the
Fejér–Riesz constraint that depends on $M[k,k] - M[K+k,K+k]$.

In v47 we add an **upper bound** on $M[k,k]$ coming from a sharp
trigonometric Cauchy–Schwarz inequality:
$$
a_k^2 \;=\; \Bigl(\int f \cos(2\pi k x)\,dx\Bigr)^2
\;\stackrel{\text{C–S}}{\le}\; \Bigl(\int f \cdot 1\Bigr)\Bigl(\int f \cos^2(2\pi k x)\Bigr)
\;=\; 1 \cdot \frac{1 + a_{2k}}{2}
\;=\; \frac{1 + a_{2k}}{2}.
$$
Similarly,
$$
b_k^2 \;\le\; \frac{1 - a_{2k}}{2},
\qquad
4\,a_k^2 b_k^2 \;\le\; (1 + a_{2k})(1 - a_{2k}) \;=\; 1 - a_{2k}^2.
$$

## New constraints

Lifting $a_k^2 \to M[k,k]$, $b_k^2 \to M[K+k,K+k]$, $a_k b_k \to M[k,K+k]$:

**(C1)** $M[k,k] \;\le\; \dfrac{1 + a_{2k}}{2}$ for $k = 1, \dots, K$.

**(C2)** $M[K+k, K+k] \;\le\; \dfrac{1 - a_{2k}}{2}$ for $k = 1, \dots, K$.

**(C3)** $4\,M[k,K+k]^2 \;\le\; 1 - a_{2k}^2$ for $k = 1, \dots, K$.
    *(SOC; or implied by C1, C2 plus Schur from $M \succeq 0$.)*

These require auxiliary variables $a_{K+1}, \dots, a_{2K}$ with their
own cell bounds (as in v46) and the extended Bochner–Toeplitz
$(2K+1)$ constraint.

## Validity

(C1)–(C3) are direct consequences of:

1. **Cauchy–Schwarz** for the inner product
   $\langle g, h\rangle_f := \int f\,g\,h$.  Since $f \ge 0$,
   this is a positive semi-definite bilinear form, and
   $\langle g, h\rangle_f^2 \le \langle g, g\rangle_f \cdot \langle h, h\rangle_f$.
2. The **trigonometric identities**
   $\cos^2 \theta = (1 + \cos 2\theta)/2$ and
   $\sin^2 \theta = (1 - \cos 2\theta)/2$ and
   $\cos\theta \sin\theta = \sin(2\theta)/2$.

Applying C–S with $g = 1, h = \cos(2\pi k x)$:
$\langle 1, \cos\rangle_f^2 = a_k^2$, $\langle 1, 1\rangle_f = 1$,
$\langle \cos, \cos\rangle_f = \int f \cos^2 = \int f (1+\cos 2kx)/2 = (1+a_{2k})/2$.
Hence $a_k^2 \le (1+a_{2k})/2$, giving (C1) at the rank-1 level.

For the lifted $M[k,k]$: we need to show that $M[k,k] \le (1+a_{2k})/2$
holds at every rank-1 feasible $M$.  At rank-1, $M[k,k] = a_k^2$, and
the C–S inequality gives $a_k^2 \le (1+a_{2k})/2$, hence
$M[k,k] \le (1+a_{2k})/2$.  $\square$

**Extension to non-rank-1 $M$:** The constraint is added as a
separate linear inequality, and is satisfied by every rank-1
feasible point.  Hence it is a **valid** relaxation constraint to add
to the SDP, and adding it can only **tighten** the relaxation (i.e.,
*increase* the optimal $\Omega$).  $\square$

The auxiliary variables $a_{2k}$ for $k = 1, \dots, K$ (i.e.,
$a_2, a_4, \dots, a_{2K}$) are needed.  We add all of $a_{K+1},
\dots, a_{2K}$ and $b_{K+1}, \dots, b_{2K}$ for symmetry, with their
cell bounds from $p$, just as in v46.

(C3) follows similarly: $a_k b_k = \langle \cos, \sin\rangle_f$ (with
inner product weighted by $f$), C–S gives
$(a_k b_k)^2 \le (1+a_{2k})/2 \cdot (1-a_{2k})/2 = (1-a_{2k}^2)/4$.
Lifting $a_k b_k \to M[k,K+k]$:
$4 M[k,K+k]^2 \le 1 - a_{2k}^2$.  This is a rotated SOC constraint
in CVXPY.

## Implementation note

(C1) and (C2) are simple linear inequalities.  (C3) is rotated SOC,
but it is **implied** by (C1) and (C2) combined with the Schur
constraint $M[k,K+k]^2 \le M[k,k] M[K+k,K+k]$ (which comes from
$M \succeq 0$):
$$
M[k,K+k]^2 \;\le\; M[k,k] M[K+k,K+k] \;\le\; \frac{(1+a_{2k})(1-a_{2k})}{4}
= \frac{1 - a_{2k}^2}{4}.
$$
So (C3) is **redundant** given (C1), (C2), and $M \succeq 0$, and
need not be added.  Only (C1) and (C2) are new.

## Construction

v47 = v45 + v46 (extended $a, b$ up to $2K$ with cell bounds and
$(2K+1)$-Bochner–Toeplitz) + the C–S constraints (C1), (C2).

## Expected gain

Hard to predict.  The C–S constraints couple $M[k,k]$ (a "free"
variable in v45) directly to $a_{2k}$ (a constrained variable via
cell bounds and Bochner).  This should significantly tighten the
moment lift, especially at moderate-to-large $K$ where $a_{2k}$ is
informative.

## Validity caveat

Note: (C1) is $M[k,k] \le (1+a_{2k})/2$.  In the SDP, $M[k,k]$ can
be DRIVEN HIGH (by the optimizer) to weaken the Fejér–Riesz
constraint.  Now it is bounded above by $(1+a_{2k})/2 \le 1$.  This
could materially change the optimum.

## Status

Awaiting theorist verdict.
