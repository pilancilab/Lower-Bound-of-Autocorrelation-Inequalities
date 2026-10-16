# v8 — Rigorous proof (and a load-bearing caveat)

## Goal

v8 attempts the verifier-recommended F-side / Fourier-side relaxation:
introduce both f-side cell averages $p_j$ and F-side cell averages $W_m$,
plus Fourier coefficients $a_k$ of $f$ with the SDP lift on $M[k,k] \ge a_k^{2}$,
and link to F via the squared identity $\hat F(k) = \hat f(k)^{2}$ giving
$A_k(F) = a_k^{2}$ for even $f$.

## Setup

Even-restricted: $f \ge 0$, $\operatorname{supp} f \subset [-1/4, 1/4]$,
$\int f = 1$, $f$ even.  $F = f * f$, $\Omega = \|F\|_\infty$.
Discretization $L_f = L_F = 1/(2N)$ on each side.

## Constraints (rigorous)

(P0) $\sum p = 1$, $L_F \sum W = 1$, $W_m \le \Omega$.
(P1) Even symmetry: $p_j = p_{2N+1-j}$, $W_m = W_{2N+1-m}$.
(P2) f-side bounds:
  $\sum_j p_j (\min\cos)_{j,k} \le a_k \le \sum_j p_j (\max\cos)_{j,k}$.
(P3) PSD lift on $M$: $M[0,0]=1$, $M[0,k]=a_k$, $M \succeq 0$.
(P5) F-side Fourier link:
  $L_F \sum_m W_m (\min\cos)_{m,k} \le M[k,k] \le L_F \sum_m W_m (\max\cos)_{m,k}$.

## NON-rigorous constraint (and why)

(P4) **Naive Parseval lower bound** (non-rigorous as written):
  $\sum_{k=1}^{K} M[k,k] \ge 1/2$.

The motivation is the Parseval identity for even $f$:
$\sum_{k=1}^{\infty} a_k^{2} = (\int f^{2} - 1)/2 \ge 1/2$ by Cauchy-Schwarz
on the support of $f$.

**Caveat (acknowledged failure of rigor).**  For finite $K$, this inequality
holds only for the **infinite** sum.  The truncation
$\sum_{k=1}^{K} a_k^{2}$ can be strictly less than $1/2$ even for admissible
$f$ — for example, the uniform $f = 2 \cdot \mathbf 1_{[-1/4, 1/4]}$ gives

$$
\sum_{k=1}^{4} a_k^{2}
\;=\; \sum_{\substack{k=1\\k\text{ odd}}}^{4} \!\!\Bigl(\frac{2}{\pi k}\Bigr)^{\!2}
\;=\; \frac{4}{\pi^{2}} + \frac{4}{9 \pi^{2}}
\;\approx\; 0.405 + 0.045 = 0.45 < 0.5,
$$

so the constraint $\sum_{k=1}^{4} M[k,k] \ge 0.5$ **excludes** the
admissible uniform $f$.  Therefore (P4) is NOT a valid relaxation
constraint at finite $K$, and v8's claim
$\Omega^{(8)} \le C^{\mathrm{even}}_{6.2}$ is **not certified** as written.

The rigorous version would require an explicit tail bound
$\sum_{k > K} a_k^{2} \le T_K$ for some computable $T_K$, but no such
bound is available because $\int f^{2}$ is unbounded above on the
admissible set (e.g., $f = c \cdot \mathbf 1_{[-1/(4c), 1/(4c)]}$ for
$c \to \infty$ has $\int f^{2} \to \infty$).

## Numerical observation

At the chosen finite $(N, K)$, MOSEK returns $\Omega^{(8)} = 1$ (the
trivial bound).  This is consistent with the relaxation finding a
non-uniform $(p, W, M)$ that satisfies the (non-rigorous) Parseval
constraint while keeping $\Omega = 1$.  Even ignoring the rigor issue,
the constraint is not binding the way I had hoped.

## Conclusion: v8 is INVALID

* The Parseval constraint as stated (P4) is invalid for finite $K$
  because it excludes admissible $f$ such as uniform $f$.
* Without (P4), v8 reduces to a strictly weaker version of v7 (since
  the F-side variables $W_m$ are loosely linked to the f-side via the
  Fourier link, which the relaxation absorbs).

v9 will:
1. Replace the naive Parseval with a **Fejér kernel test functional**
   that gives a rigorous linear inequality without any tail issue.
2. Use the cleaner $\int F_K(t) F(t)\,dt \le \Omega \cdot \int F_K(t)\,dt = \Omega$
   from the verifier's recommendation, where $F_K$ is the (positive)
   Fejér kernel of order $K$.
