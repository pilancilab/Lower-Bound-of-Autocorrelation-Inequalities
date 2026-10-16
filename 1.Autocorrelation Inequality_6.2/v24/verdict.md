# Verdict — Autocorrelation Inequality 6.2, v24

## 1. Summary

v24 introduces a new LP variable \(M_\infty \ge 0\) representing an
upper bound on \(\|f\|_\infty\), with three constraints:

* **(M1)** \(M_\infty \ge p_j/L_f\) for each \(j\).
* **(M2)** \(M_\infty \ge \Omega\).
* **(M3)** \(M_\infty \ge 2\).

All three are individually valid, but the LP optimum is identical
to v14 (no improvement). The proposer correctly identifies the
vacuousness.

## 2. Verdict

**VALID** (rigorous but vacuous).

### Verification

* **(M1)** \(p_j/L_f\) is the cell-average density of \(f\) on
  \(I_j\), bounded by \(\|f\|_\infty\). So \(M_\infty\ge\|f\|_\infty\ge p_j/L_f\). ✓
* **(M2)** \(\Omega = \|f*f\|_\infty \le \|f\|_2^{\,2} \le \|f\|_\infty\cdot\|f\|_1 = \|f\|_\infty\),
  using Young's \(\|h*h\|_\infty\le\|h\|_2^{\,2}\) and Cauchy-Schwarz.
  So \(M_\infty\ge\|f\|_\infty\ge\Omega\). ✓
* **(M3)** \(\|f\|_\infty\ge \int f/|\operatorname{supp}f| = 1/(1/2) = 2\)
  for \(f\ge 0\) supported in \([-1/4,1/4]\) with \(\int f = 1\). So
  \(M_\infty\ge 2\). ✓ (The proof's "Hölder" derivation in lines
  24–28 is muddled, but the parenthetical "for uniform \(f = 2\),
  \(\|f\|_\infty = 2\), and any concentration only increases
  \(\|f\|_\infty\)" gives the correct argument.)

### Vacuousness

The LP minimizes \(\Omega\) subject to constraints. \(M_\infty\)
is a free variable that doesn't appear in the objective. The LP
can set \(M_\infty\) to any value \(\ge \max(p_j/L_f, \Omega, 2)\),
which is satisfiable for any feasible \((p, \Omega)\) with \(\Omega\le M_\infty\).

So **\(M_\infty\) does not constrain the LP at all**: the LP picks
\((p, \Omega)\) to minimize \(\Omega\) subject to the v14 constraints
(Fejér + cell-bounds), then sets \(M_\infty\) to whatever satisfies
the new constraints. Adding (M1), (M2), (M3) doesn't change the
minimum of \(\Omega\).

## 3. Comments

### Comment 1. Same vacuousness pattern as v15, v16, v18, v20.

| iteration | new content | result |
|---|---|---|
| v15 | shifted (W) | no change |
| v16 | Bochner Toeplitz | no change |
| v18 | L³ constraint | no change |
| v19 | multi-\(K\) Fejér | tiny change |
| v20 | Lasserre level-2 | no change |
| v21 | squared Fejér kernel | tiny change |
| **v24** | **\(M_\infty\) upper-bound variable** | **no change** |

The pattern continues: any constraint that doesn't directly enter
the objective \(\Omega\) (or doesn't tighten an existing binding
constraint) is structurally vacuous in the v14 minimization.

### Comment 2. The proposer correctly diagnoses this.

The proposer's conclusion:
> $M_\infty$ records the f-side $L^\infty$ bound but is not in any
> binding constraint that pushes $\Omega$ up.
> 
> The persistent observation across v15–v24: adding more constraints
> to the v14 framework does not improve the bound at small $N, K$.
> The bottleneck is fundamentally the rigorous cell-bound on $a_k$
> from $p$, which is $O(L_f \cdot k)$ and slow to converge as $N$
> grows.

This is the right diagnosis. The path forward (per v20) is
parameter scaling or a structurally different construction, not
more auxiliary variables.

---

**Final verdict: VALID** (no improvement over v14).
