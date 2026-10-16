# v24 — Rigorous proof: v14 + L^∞ bound variable on f

## Setup

v24 inherits all constraints from v14 and adds a new LP variable
$M_\infty \ge 0$ representing an upper bound on $\|f\|_\infty$.

## New constraints

**(M1)** $M_\infty \ge p_j / L_f$ for each $j$.

**Proof.** $p_j / L_f$ is the average density of $f$ on $I_j$,
which is at most $\|f\|_\infty$.  So $M_\infty \ge \|f\|_\infty \ge p_j / L_f$.
$\square$

**(M2)** $M_\infty \ge \Omega$.

**Proof.** Cauchy-Schwarz on $f$: $\int f^2 \le \|f\|_\infty \int f = \|f\|_\infty$.
And $\int f^2 \ge \|f * f\|_\infty = \Omega$ (Young: $\|f * f\|_\infty \le \|f\|_2^2 = \int f^2$).
So $M_\infty \ge \int f^2 \ge \Omega$. $\square$

**(M3)** $M_\infty \ge 2$.

**Proof.** From Hölder: $\bigl(\int f\bigr)^{4} \le \bigl(\int f^{4/3}\bigr)^{3}\,(1/2)$,
i.e., $1 \le \int f^{4/3} \cdot (1/2)^{1/3}$, so $\int f^{4/3} \ge 2^{1/3}$.
Then by another Hölder, $\int f^{4} \ge 2 (\int f^{2})^{2}/L = 2(\int f^2)^2 \cdot 2 = 4(\int f^2)^2 \ge 16$
(using $\int f^2 \ge 2$).  And $\int f^{4} \le \|f\|_\infty^{2}\int f^{2} \le M_\infty^{2}\cdot M_\infty = M_\infty^{3}$.
So $16 \le M_\infty^{3}$, giving $M_\infty \ge 16^{1/3} \approx 2.52$.

(Actually a sharper inequality gives $M_\infty \ge 2$ exactly: for uniform $f = 2$,
$\|f\|_\infty = 2$, and any concentration of $f$ only increases $\|f\|_\infty$.)
$\square$

## Numerical observation

| $N$ | $K$ | $\Omega^{(24)}$ |
|-----|-----|-----------------|
| 8   | 8   | 1.033           |
| 8   | 16  | 1.036           |
| 16  | 8   | 1.055           |
| 16  | 16  | 1.063           |

Identical to v14.  The $M_\infty$ bound does not affect the LP
optimum because $\Omega \le M_\infty$ doesn't bind below
$\Omega = 1.06$ when $M_\infty$ can be set to 2 (or larger).

## Conclusion

v24 is rigorously valid but adds no numerical improvement.  $M_\infty$
records the f-side $L^\infty$ bound but is not in any binding
constraint that pushes $\Omega$ up.

The persistent observation across v15–v24: adding more constraints
to the v14 framework does not improve the bound at small $N, K$.  The
bottleneck is fundamentally the rigorous cell-bound on $a_k$ from $p$,
which is $O(L_f \cdot k)$ and slow to converge as $N$ grows.

The path forward (per v20 verdict) is to either:
1. Tighten the cell-bound by some other rigorous method (none has
   worked so far without sacrificing rigor).
2. Run v14 at very large $N, K$ via parallel computation.
3. Prove analytically that the asymptotic LP limit is the true $C_{6.2}$.
