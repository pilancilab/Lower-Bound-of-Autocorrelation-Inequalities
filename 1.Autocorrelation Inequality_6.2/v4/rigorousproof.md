# v4 — Rigorous proof of the multi-window constraints

## Setup

We continue with the v0/v1 (p, Q) SDP lift and add a new family of
constraints: **central-window integral inequalities at multiple half-widths**
plus **shifted-window inequalities**.

Notation: $f \ge 0$, $\operatorname{supp} f \subset [-1/4, 1/4]$,
$\int f = 1$, $F = f * f$, $\Omega = \|F\|_\infty$.  Discretization:
$2N$ intervals $I_j = [-1/4 + (j-1)L,\, -1/4 + jL]$ of width $L = 1/(4N)$.
Variables: $p_j = \int_{I_j} f$, $Q[j,k] = p_j p_k$.

## New constraint family (W) — windowed integral bounds

For any window $B = [a,\, b] \subset [-1/2, 1/2]$ with width $h_B = b - a$,

$$
\int_B (f*f)(t)\,dt
\;=\; \iint f(\tau)\,f(x)\,\mathbf 1_{a \le \tau + x \le b}\,d\tau\,dx
\;\le\; h_B \cdot \Omega.
\qquad (\dagger)
$$

The first equality uses Fubini and the substitution $x = t - \tau$ as in
the v1 (C) proof.  The inequality uses $F \le \Omega$.

**Claim.**  For every admissible $f$,

$$
\sum_{(j,k):\, I_j \times I_k \subset \{(\tau,x):\, a \le \tau + x \le b\}}\!\!\!\! Q[j,k]
\;\le\; h_B \cdot \Omega.
$$

**Proof.**  By the same argument as in v1 (C):
$$
\iint f(\tau) f(x)\,\mathbf 1_{a \le \tau+x \le b}\,d\tau\,dx
\;\ge\; \sum_{(j,k):\, I_j \times I_k \subset \text{preimage}} \,p_j p_k
\;=\; \sum_{(j,k):\, I_j \times I_k \subset \text{preimage}} Q[j,k],
$$
because on each fully-contained cell the indicator is identically $1$ and
$\iint_{I_j \times I_k} f(\tau)f(x)\,d\tau\,dx = p_j p_k$ exactly.
Combining with $(\dagger)$ proves the claim. $\square$

The cell $I_j \times I_k$ is fully contained in
$\{(\tau,x): a \le \tau+x \le b\}$ iff
$$
a \;\le\; -1/2 + (j+k-2)L \quad\text{and}\quad -1/2 + (j+k)L \;\le\; b,
$$
i.e.,
$$
\frac{2(a+1/2)}{L} + 2 \;\le\; j+k \;\le\; \frac{2(b+1/2)}{L}.
$$
We add a windowed constraint for each $B$ in a chosen family.

## Choice of windows in v4

We add constraints for the following windows on top of the v1 (C)
single central window:

* Centered windows $B = [-h, h]$ for $h$ varying through several
  half-widths, namely $h \in \{2L, 3L, \dots, 1/2\}$.
* Shifted windows $B = [c - h, c + h]$ for $c \in \{-1/4, +1/4\}$ and
  the same set of half-widths.

Each windowed constraint is a linear inequality on $Q$ derived above and
is therefore **valid for every admissible $f$**.

## Validity of the relaxation

For every admissible $f$ the rank-one $Q^* = pp^{\!\top}$ satisfies the
v0/v1 constraints (PSD lift, mass, marginals, $Q \ge 0$, bands,
central-window) and every windowed constraint (W) by the claim above.
Therefore the v4 optimum
$$
\Omega^{(4)} \;\le\; \|f * f\|_\infty
$$
for every admissible $f$, giving $\Omega^{(4)} \le C_{6.2}$. $\blacksquare$

## Numerical observation

At $N = 8$, MOSEK still returns $\Omega^{(4)} = 1$.  The reason is the
same diagonal-$Q$ degeneracy that defeated v1, v2, v3: the matrix
$Q = (1/(2N))\,I$ together with $p = (1/(2N))\,\mathbf 1$ satisfies every
linear and PSD constraint above, including all windowed constraints,
because:

* For a window $B$ containing $D_B$ diagonal cells $(j,j)$, the diagonal
  $Q$ contributes $D_B / (2N)$ to the windowed sum.
* Counting: $D_B$ scales linearly with the window length $h_B$ as
  $D_B \approx 2 h_B \cdot N$, so
  $\sum_{D_B} Q[j,j] = 2 h_B \cdot N \cdot 1/(2N) = h_B$.
* Hence the constraint $\sum \le 2 h_B \cdot \Omega$ becomes
  $h_B \le 2 h_B \cdot \Omega$, i.e., $\Omega \ge 1/2$ — strictly weaker
  than the trivial $\Omega \ge 1$ from the band sums.

So the windowed constraints, while rigorous, are **strictly redundant**
for the diagonal $Q$.  The diagonal solution evades them entirely and
gives $\Omega = 1$.

This confirms the diagnosis from v1 and v3: linear and PSD relaxations
on $(p, Q)$ alone cannot exclude the basis-vector mixture, and the
relaxation cannot strictly improve on $\Omega = 1$ without introducing
constraints that fundamentally distinguish a single $f$ from a mixture of
delta-supported $f$'s.

v5 will attempt this distinction by working **on the dual side**: rather
than relax the primal, we exhibit a dual feasible point (a positive
measure $\nu$ on $[-1/2, 1/2]$) for which
$\inf_f \int (f*f) d\nu / \int d\nu > 1$, giving an honest lower bound.
