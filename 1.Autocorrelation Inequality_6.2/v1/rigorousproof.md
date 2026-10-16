# v1 — Rigorous proof of new constraints

We work in the normalized setting

> $f \ge 0$,  $\operatorname{supp}(f) \subset [-1/4, 1/4]$, $\int f = 1$,
> $F(t) = (f*f)(t)$,  $\Omega = \|F\|_\infty$,

and minimize $\Omega$.  The discretization uses $2N$ intervals
$I_j = [-1/4 + (j-1)L,\; -1/4 + jL]$, $j = 1,\dots,2N$, of equal width
$L = 1/(4N)$.  The lifting variables are

$$
p_j \;:=\; \int_{I_j} f \;\ge\; 0, \qquad
Q[j,k] \;:=\; p_j p_k \quad (1 \le j,k \le 2N).
$$

In v0 the only constraints on $Q$ were
$\bigl[\begin{smallmatrix}1 & p^{\!\top}\\ p & Q\end{smallmatrix}\bigr] \succeq 0$
(i.e. $Q \succeq pp^{\!\top}$ in Loewner order),
plus the band-integral inequalities derived from $F \le \Omega$.  v1 adds
three constraints, each of which holds for every true $Q = pp^{\!\top}$
coming from an admissible $f$.

---

## Constraint (A) — element-wise nonnegativity

**Claim.**  $Q[j,k] \ge 0$ for all $1 \le j,k \le 2N$.

**Proof.**  $f \ge 0$ implies $p_j = \int_{I_j} f \ge 0$ for every $j$.  Hence
$Q[j,k] = p_j p_k \ge 0$ as a product of two nonnegative numbers. $\square$

---

## Constraint (B) — row-sum (marginal) constraint

**Claim.**  $\sum_{k=1}^{2N} Q[j,k] = p_j$ for all $j$.

**Proof.**  By definition $Q[j,k] = p_j p_k$ and the program already enforces
$\sum_{k=1}^{2N} p_k = 1$ (the normalization $\int_{-1/4}^{1/4} f = 1$
discretized).  Hence
$$
\sum_{k=1}^{2N} Q[j,k] \;=\; p_j \sum_{k=1}^{2N} p_k \;=\; p_j \cdot 1 \;=\; p_j. \qquad\square
$$

---

## Constraint (C) — central-window integral constraint

This is the **key new constraint** of v1, the analogue of integrating
$F = f*f$ against the indicator of $[-1/4, 1/4]$ on the $t$-axis.

**Claim.**  Let $\mathcal C$ denote the set of index pairs $(j,k)$ such that
the cell $I_j \times I_k$ is entirely contained in
$\{(\tau, x) \in [-1/4,1/4]^2 \;:\; |\tau + x| \le 1/4\}$.  Equivalently,

$$
\mathcal C \;=\; \{(j,k) \;:\; N+2 \le j+k \le 3N\}.
$$

Then for every admissible $f$,

$$
\sum_{(j,k)\in\mathcal C} Q[j,k] \;\le\; \frac{\Omega}{2}.
$$

**Proof.**  Two facts.  First, since $F \ge 0$ and $F \le \Omega$ on
$[-1/2, 1/2]$,
$$
\int_{-1/4}^{1/4} F(t)\,dt \;\le\; \Omega \cdot \int_{-1/4}^{1/4} dt \;=\; \frac{\Omega}{2}.
\tag{1}
$$
Second, applying Fubini to the convolution and then Tonelli to the indicator,
$$
\int_{-1/4}^{1/4} F(t)\,dt
\;=\; \int_{-1/4}^{1/4}\!\!\int_{\mathbb R}\! f(\tau)\,f(t-\tau)\,d\tau\,dt
\;=\; \iint f(\tau)f(x)\,\mathbf 1_{|\tau+x|\le 1/4}\,d\tau\,dx,
$$
where the substitution $x = t - \tau$ exchanged the order of integration.

The double integral splits over the cells $I_j \times I_k$:
$$
\iint f(\tau)f(x)\,\mathbf 1_{|\tau+x|\le 1/4}\,d\tau\,dx
\;=\; \sum_{j,k} \iint_{I_j\times I_k} f(\tau)f(x)\,\mathbf 1_{|\tau+x|\le 1/4}\,d\tau\,dx.
\tag{2}
$$

For each cell index $(j,k)$ we estimate the inner integral.  Let
$\chi_{j,k}(\tau,x) := \mathbf 1_{|\tau+x|\le 1/4}$ restricted to
$(\tau,x) \in I_j \times I_k$.  By construction,
$\chi_{j,k} \equiv 1$ on the entire cell exactly when
$(j,k) \in \mathcal C$, and $\chi_{j,k} \ge 0$ everywhere.  Therefore
$$
\iint_{I_j\times I_k} f(\tau)f(x)\,\chi_{j,k}(\tau,x)\,d\tau\,dx
\;\ge\; \mathbf 1_{(j,k)\in\mathcal C}\, \iint_{I_j\times I_k} f(\tau)f(x)\,d\tau\,dx,
$$
because we may drop the indicator $\chi_{j,k}$ outright on cells in
$\mathcal C$ and bound the integral by zero (the trivial lower bound)
otherwise.  By Fubini on the rectangle $I_j \times I_k$,
$$
\iint_{I_j\times I_k} f(\tau)f(x)\,d\tau\,dx
\;=\; \Bigl(\int_{I_j} f\Bigr)\Bigl(\int_{I_k} f\Bigr) \;=\; p_j p_k \;=\; Q[j,k].
$$
Substituting back into (2),
$$
\int_{-1/4}^{1/4} F(t)\,dt
\;\ge\; \sum_{(j,k)\in\mathcal C} Q[j,k].
\tag{3}
$$

Combining (1) and (3) gives the claim:
$\sum_{\mathcal C} Q[j,k] \le \Omega/2$. $\square$

**Identification of $\mathcal C$.**  A cell $I_j \times I_k$ is contained in
$\{|\tau+x| \le 1/4\}$ iff
$$
\min_{(\tau,x)\in I_j\times I_k} (\tau+x) \;=\; -1/2 + (j+k-2)L \;\ge\; -1/4
\quad\text{and}\quad
\max_{(\tau,x)\in I_j\times I_k} (\tau+x) \;=\; -1/2 + (j+k)L \;\le\; 1/4.
$$
Substituting $L = 1/(4N)$ gives $(j+k-2)/(4N) \ge 1/4$ and $(j+k)/(4N) \le 3/4$,
i.e. $j+k \ge N+2$ and $j+k \le 3N$, as claimed.

---

## Validity of the relaxation

For any admissible $f$, the assignment
$$
\Omega^* \;=\; \|f * f\|_\infty,
\quad p_j^* \;=\; \int_{I_j} f,
\quad Q^*[j,k] \;=\; p_j^* p_k^*
$$
is feasible for the v1 program.  Indeed:

* $p^* \ge 0$ and $\sum p^* = 1$ are immediate.
* The PSD-lift $\bigl[\begin{smallmatrix}1 & p^{*\top}\\ p^* & Q^*\end{smallmatrix}\bigr]$
  is rank-one PSD.
* Constraints (A), (B), (C) hold by the three claims above.
* The v0 band-integral constraints are derived in the same way as (1)–(3)
  but with the indicator of $B_m = [-1/2 + (m-1)L,\, -1/2 + mL]$ instead of
  $[-1/4, 1/4]$, and they are still valid.

Hence the v1 optimum $\Omega^{(1)} \le \Omega^*$.  Since this holds for
*every* admissible $f$,
$$
\Omega^{(1)} \;\le\; \inf_{f}\, \|f * f\|_\infty \;=\; C_{6.2},
$$
so $\Omega^{(1)}$ is a valid lower bound on $C_{6.2}$. $\blacksquare$

---

## Numerical observation (informational, not part of the proof)

After implementing constraints (A), (B), (C) on top of v0, a MOSEK solve at
$N = 8$ still yields $\Omega^{(1)} = 1$.  The reason is that the relaxation
admits a degenerate feasible point of the form
$Q = (1/(2N))\,\mathrm I$ together with $p = (1/(2N))\,\mathbf 1$, on which
the band constraints become tight (the alternating-anti-diagonal pattern
$a_s = 1/(2N) \cdot \delta_{s\,\text{even}}$) and constraint (C) is exactly
$1/2 \le \Omega/2$ at $\Omega = 1$.  This $Q$ does not correspond to any
rank-one outer product but is still admitted by $Q \succeq pp^{\!\top}$
because $Q - pp^{\!\top} = (1/(2N))\,\mathrm I - (1/(2N))^2 \mathbf 1\mathbf 1^{\!\top}$
is PSD.  Subsequent versions will add constraints that exclude this
diagonal-only relaxation point — for instance constraints that
distinguish $Q$ from $\operatorname{diag}(p)$.

The proof of the validity of (A), (B), (C) remains correct regardless: every
true rank-one $Q^* = pp^{\!\top}$ continues to be feasible, and so the
relaxed optimum is still a valid (though not yet improved) lower bound on
$C_{6.2}$.
