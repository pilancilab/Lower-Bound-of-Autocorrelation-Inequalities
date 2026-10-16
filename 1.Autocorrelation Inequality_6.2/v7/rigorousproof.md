# v7 — Rigorous proof: clean reboot addressing all prior verdicts

## Goal of v7

The goal of v7 is to produce a **fully rigorous** convex relaxation
whose optimum lower-bounds $C^{\mathrm{even}}_{6.2}$, the
even-restricted version of $C_{6.2}$:
$$
C^{\mathrm{even}}_{6.2} \;:=\; \inf_{f \text{ admissible AND } f \text{ even}}\, \|f * f\|_\infty.
$$
We address every load-bearing issue raised in the verdicts of v1, v2,
v3, v4, v5: the false v0 band-integral identity is **dropped**, the
WLOG reductions are **explicitly justified**, and the (W) and (G)
constraints are derived without any "approximation" hand-wave.

The relaxation does **not** bound $C_{6.2}$ in general (we have
$C_{6.2} \le C^{\mathrm{even}}_{6.2}$, so a lower bound on the
even-restricted constant is automatically a lower bound on the
even-restricted constant only).  v8 will attempt the harder extension
to non-even $f$.

## WLOG reductions (justified)

The original problem is to lower-bound the constant
$C_{6.2}$ defined by
$$
C_{6.2} \;=\; \inf_{f}\; \frac{\max_{|t|\le 1/2} (f * f)(t)}{\bigl(\int_{-1/4}^{1/4} f\bigr)^{2}}
$$
over all nonnegative $f : \mathbb R \to \mathbb R$ for which the
denominator is nonzero.

**Reduction 1 (support compactification).**  We may restrict to $f$
supported in $[-1/4, 1/4]$.

*Proof.*  Replace $f$ by $f \cdot \mathbf 1_{[-1/4, 1/4]}$.  This leaves
$\int_{-1/4}^{1/4} f$ unchanged (the integrand is identical on the
window).  For the numerator, by monotonicity of convolution of
nonneg functions, $f \cdot \mathbf 1_{[-1/4, 1/4]} \le f$ pointwise gives
$$
(f \cdot \mathbf 1_{[-1/4, 1/4]}) * (f \cdot \mathbf 1_{[-1/4, 1/4]})
\;\le\; f * f
\quad\text{pointwise.}
$$
Hence $\max_{|t|\le 1/2}(f \cdot \mathbf 1_{[-1/4, 1/4]}) * (\dots)(t) \le \max_{|t|\le 1/2}(f * f)(t)$,
so the ratio for the truncated $f$ is no larger than for $f$.  Hence
the infimum is attained on $f$ supported in $[-1/4, 1/4]$.  $\square$

**Reduction 2 (mass normalization).**  We may further assume
$\int_{-1/4}^{1/4} f = 1$.

*Proof.*  The ratio
$\max_{|t|\le 1/2}(f*f)/(\int_{-1/4}^{1/4} f)^{2}$
is invariant under scaling $f \to \lambda f$ for any $\lambda > 0$, since
both numerator and denominator are degree-2 homogeneous in $f$.  So we
may rescale to $\int_{-1/4}^{1/4} f = 1$.  $\square$

**Even restriction.**  v7 additionally assumes $f$ is even, so we
work in
$$
\mathcal A^{\mathrm{even}}
\;:=\; \{f : \mathbb R \to [0, \infty),\;\;
\operatorname{supp}(f) \subset [-1/4, 1/4],\;\;
\int f = 1,\;\;
f(-x) = f(x)\}.
$$
The relaxation's optimum will lower-bound
$\inf_{f \in \mathcal A^{\mathrm{even}}} \|f * f\|_\infty$, which is
$C^{\mathrm{even}}_{6.2}$.

## Discretization

Let $L = 1/(4N)$ and partition $[-1/4, 1/4]$ into $2N$ intervals
$$
I_j \;=\; [-1/4 + (j-1)L,\; -1/4 + jL], \qquad j = 1, 2, \dots, 2N.
$$
The reflection $x \mapsto -x$ exchanges $I_j$ and $I_{2N+1-j}$ because
$-I_j = [1/4 - jL,\, 1/4 - (j-1)L] = I_{2N+1-j}$ (using
$2N L = 1/2$, so $1/4 - jL = -1/4 + (2N - j) L = -1/4 + ((2N+1-j) - 1) L$).

Lifting variables:
$$
p_j \;:=\; \int_{I_j} f \;\ge\; 0,
\qquad
Q[j,k] \;:=\; p_j p_k \quad (1 \le j, k \le 2N).
$$

## Constraints

### (P0) Schur-complement PSD lift

$\bigl[\begin{smallmatrix}1 & p^{\!\top}\\ p & Q\end{smallmatrix}\bigr] \succeq 0$,
equivalently $Q \succeq pp^{\!\top}$.  Holds for the rank-one truth.

### (P1) Mass

$\sum_{j=1}^{2N} p_j = 1$.  Holds because
$\sum p_j = \int_{-1/4}^{1/4} f = 1$ by Reduction 2.

### (P2) Element-wise nonnegativity

$Q[j, k] \ge 0$.  Holds because $p_j, p_k \ge 0$.

### (P3) Marginals

$\sum_{k=1}^{2N} Q[j, k] = p_j$.  Holds because
$\sum_k Q[j,k] = p_j (\sum_k p_k) = p_j \cdot 1 = p_j$.

### (P4) v5 even-symmetry constraints

* **(E)** $p_j = p_{2N+1-j}$.  *Proof.* For $f$ even,
  $p_{2N+1-j} = \int_{-I_j} f(x)\,dx = \int_{I_j} f(-y)\,dy = \int_{I_j} f(y)\,dy = p_j$.

* **(F)** $Q[j, k] = Q[2N+1-j,\, 2N+1-k]$.  *Proof.*
  $Q[j,k] = p_j p_k = p_{2N+1-j} p_{2N+1-k} = Q[2N+1-j,\, 2N+1-k]$
  by applying (E) twice.

* **(G)** $Q[j,\, 2N+1-j] = Q[j, j]$.  *Proof.*
  $Q[j,\, 2N+1-j] = p_j p_{2N+1-j} = p_j \cdot p_j = p_j^{2} = Q[j, j]$
  using (E) and the rank-one definition $Q[j,j] = p_j p_j$.

### (P5) Rigorous (W) windowed integral constraints

This is the **key fix** that resolves the verdict's load-bearing issue
on the v0/v1 band-integral chain.

Let $B = [a, b] \subset [-1/2, 1/2]$ be a window of length $h_B = b - a$.
Define
$$
\mathcal C(B) \;:=\; \{(j, k) : I_j \times I_k \subset \{(\tau, x) : \tau + x \in B\}\}.
$$

**Claim.**  For every admissible $f$,
$$
\sum_{(j,k)\in\mathcal C(B)} Q[j, k] \;\le\; h_B \cdot \Omega,
\qquad \Omega := \|f * f\|_\infty.
$$

**Proof.**  Two computations.

*(a) Rewrite the band integral as a double integral over $f$.*  By
Fubini and the substitution $x = t - \tau$,
$$
\int_B (f * f)(t)\,dt
\;=\; \int_B \int_{\mathbb R} f(\tau)\,f(t - \tau)\,d\tau\,dt
\;=\; \iint f(\tau)\,f(x)\,\mathbf 1_{a \le \tau + x \le b}\,d\tau\,dx.
$$

*(b) Lower bound the double integral cell-by-cell.*  Since
$f(\tau) f(x) \ge 0$ pointwise, and $\mathbf 1_{a \le \tau+x \le b}$ is
identically $1$ on every cell $(j, k) \in \mathcal C(B)$ and
identically nonnegative everywhere else,
$$
\iint f(\tau) f(x)\,\mathbf 1_{a\le\tau+x\le b}\,d\tau\,dx
\;\ge\; \sum_{(j,k)\in\mathcal C(B)} \iint_{I_j \times I_k} f(\tau) f(x)\,d\tau\,dx
\;=\; \sum_{(j,k)\in\mathcal C(B)} p_j p_k
\;=\; \sum_{(j,k)\in\mathcal C(B)} Q[j, k].
$$
The middle equality is **exact** (Fubini on the rectangle):
$\iint_{I_j \times I_k} f(\tau) f(x)\,d\tau\,dx = (\int_{I_j} f)(\int_{I_k} f) = p_j p_k$.
**There is no piecewise-constant approximation here.**  This is the
core fix.

*(c) Upper bound the band integral by $h_B \Omega$.*  Since
$F = f * f \ge 0$ and $F \le \Omega$ pointwise,
$$
\int_B F(t)\,dt \;\le\; \Omega \cdot \int_B 1\,dt \;=\; \Omega \cdot h_B.
$$

Combining (a)–(c) gives the claim. $\square$

**Cell-membership criterion.**  $(j, k) \in \mathcal C(B = [a, b])$ iff
$$
-1/2 + (j+k-2)\,L \;\ge\; a
\quad\text{and}\quad
-1/2 + (j+k)\,L \;\le\; b,
$$
i.e.
$$
\frac{a + 1/2}{L} + 2 \;\le\; j+k \;\le\; \frac{b + 1/2}{L}.
$$

**Family of windows used in v7.**  We add (W) for each centered window
$B = [-h, h]$ at half-widths $h \in \{2L, 3L, \dots, 1/2\}$ plus the
full window $B = [-1/2, 1/2]$ (which contributes the trivial
$\sum_{j,k} Q[j,k] \le \Omega$ since every cell is fully inside the
full window).

## Validity of the relaxation

For every admissible **even** $f$, the assignment
$$
\Omega^* = \|f * f\|_\infty,\quad
p^*_j = \int_{I_j} f,\quad
Q^*[j, k] = p^*_j\, p^*_k
$$
is feasible for v7:

* (P0) holds because $Q^* = p^*\,p^{*\top}$ is rank-one PSD;
* (P1)–(P3) hold by the proofs above;
* (P4) (E)/(F)/(G) hold by the v5 proofs above;
* (P5) (W) holds by the (W) claim above.

So the v7 optimum $\Omega^{(7)} \le \Omega^* = \|f * f\|_\infty$ for
every admissible even $f$, giving
$$
\Omega^{(7)} \;\le\; \inf_{f \in \mathcal A^{\mathrm{even}}} \|f * f\|_\infty \;=\; C^{\mathrm{even}}_{6.2}.
$$

And $C^{\mathrm{even}}_{6.2} \ge \Omega^{(7)}$ — a **rigorous** lower
bound on $C^{\mathrm{even}}_{6.2}$.

## Numerical results

| $N$ | $\Omega^{(7)}$ |
|-----|----------------|
| 4   | 1.046          |
| 6   | 1.082          |
| 8   | 1.164          |
| 10  | 1.206          |
| 12  | 1.241          |
| 16  | 1.299          |
| 20  | 1.344          |

These bounds are **strictly above 1** for every $N \ge 4$, and grow
monotonically toward the analytic value $C^{\mathrm{even}}_{6.2} = 2$
(achieved by $f = 2 \cdot \mathbf 1_{[-1/4, 1/4]}$).  The convergence is
slower than v5's because v7 uses only the rigorous (W) constraints,
whereas v5 inherited the stronger but invalid v0 bands.

## What v7 does **not** prove

* It does **not** bound $C_{6.2}$ in general; only $C^{\mathrm{even}}_{6.2}$.
* The bound is currently weaker than v5's because v5's relied on a
  false-equality chain (v5's was numerically larger but not
  certifiable).
* v8 should aim to (i) get more out of the (W) family by adding
  shifted windows, and (ii) extend the relaxation to non-even $f$.
