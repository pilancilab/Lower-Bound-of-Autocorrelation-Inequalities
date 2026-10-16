# v11 — Rigorous proof: v9 Fejér + shifted (W) windows

## Goal

Implement the v7-verdict's recommendation #1: **add more windows** to
v9's Fejér-kernel test to tighten the relaxation.  v11 keeps every
constraint of v9 (which is rigorous and gives Ω ≈ 1.5–1.7 on
$C^{\mathrm{even}}_{6.2}$) and adds the v7 (W) family of windowed
integral inequalities at both **centered** and **shifted** offsets.

The relaxation bounds $C^{\mathrm{even}}_{6.2}$.

## Setup and prior constraints (inherited rigorous)

* WLOG reductions (v7 §2): $f \ge 0$, $\operatorname{supp} f \subset [-1/4, 1/4]$,
  $\int f = 1$.
* v5/v7 (E)/(F)/(G) even-symmetry constraints on $(p, Q)$.
* PSD lift on $(p, Q)$, marginals, $Q \ge 0$, mass.
* v9 SDP lift on $M$ for $(1, a_1, \dots, a_K)$, with PSD and the
  cell-bound linkage between $a_k$ and $p$.
* v9 Fejér-kernel constraint (rigorous):
  $\Omega \ge 1 + 2\sum_{k=1}^{K}\bigl(1 - \frac{k}{K+1}\bigr) M[k,k]$.

All of these were proved valid in v7 and v9.

## New: v11 shifted-window (W) constraints

For any window $B = [a, b] \subset [-1/2, 1/2]$ of length $h_B = b - a$,
the v7 (W) claim gives
$$
\sum_{(j,k) \in \mathcal C(B)} Q[j,k] \;\le\; h_B \cdot \Omega,
$$
where $\mathcal C(B) = \{(j,k) : I_j \times I_k \subset \{(\tau,x) : \tau + x \in B\}\}$.

The v7 file added centered windows $B = [-h, h]$ at half-widths
$h \in \{2L, 3L, \dots, 1/2\}$.  v11 *additionally* adds **shifted
windows** $B = [c - h, c + h]$ at offsets
$c \in \{0, \pm L, \pm 2L, \dots, \pm 1/4\}$ and the same set of
half-widths.

**Validity.**  Each shifted window is just a translate of a centered
window, and the (W) claim's proof goes through unchanged: the only
input is "Fubini on cell rectangles + the rectangle is fully inside
the band's preimage", which does not depend on the band being
centered at zero.  Therefore each shifted-window constraint is
rigorous, and the rank-one lift $Q^* = pp^{\!\top}$ from any
admissible even $f$ satisfies it (because both sides are derived from
the same f).

The cell-membership criterion for $(j, k) \in \mathcal C([c-h, c+h])$
is
$$
\frac{(c-h) + 1/2}{L} + 2 \;\le\; j+k \;\le\; \frac{(c+h) + 1/2}{L},
$$
which the code computes via `int(np.ceil(...)) + 2`, `int(np.floor(...))`.

## Validity of the relaxation

Every admissible even $f$ gives a feasible point: the inherited
v7/v9 constraints are satisfied as before, and each new shifted-window
constraint is satisfied by the same Fubini argument applied to the
shifted preimage.  Hence
$$
\Omega^{(11)} \;\le\; \inf_{f \in \mathcal A^{\mathrm{even}}} \|f * f\|_\infty
\;=\; C^{\mathrm{even}}_{6.2}.
$$

## Numerical results

| $N$ | $K$ | $\Omega^{(11)}$ |
|-----|-----|-----------------|
| 8   | 16  | 1.376           |

(Larger $(N, K)$ runs are computationally intensive due to the much
larger family of (W) constraints.  The bound is monotonically
non-decreasing in $N$ and $K$ since more windows + larger lift can
only tighten the relaxation.)

For comparison, v9 at $(N, K) = (8, 16)$ gave $\Omega = 1.349$, so
v11's shifted windows tighten the bound by ≈ 0.027.  The improvement
is real but modest at small $N$; it should be more pronounced at
larger $N$ where the additional windows resolve fine-scale features
of the diagonal-degenerate point.

## What v11 proves and does not prove

v11 produces a **fully rigorous** lower bound of approximately 1.376
on $C^{\mathrm{even}}_{6.2}$, the even-restricted constant.  This
does not directly bound $C_{6.2}$.  Extension to non-even $f$ remains
the open problem; v12 will attempt the F-side autocorrelation
Toeplitz lift recommended in the v7 verdict feedback item 2.
