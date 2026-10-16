# v16 — Rigorous proof: combination of v14 + v13 + v11

## Setup

v16 combines all the rigorously valid constraints discovered in
prior iterations:

* (p, Q) PSD lift, mass, marginals, $Q \ge 0$ (from v0/v7).
* v7/v11 (W) windowed integral inequalities at multiple
  centered + shifted windows.
* (a, b) f-side Fourier coefficients with v2/v7 cell bounds from $p$.
* M PSD lift on $(1, a_1, \dots, a_K, b_1, \dots, b_K)$.
* v14 autocorrelation Fejér constraint via the $u_k = M[k,k]+M[K+k,K+k]$,
  $v_k \ge u_k^{2}$ chain.
* v13 Bochner Toeplitz PSD constraint on the f-side Fourier sequence
  $(\hat f(k))_{k=0}^{K}$.

Each constraint is individually valid for any admissible $f$, so
the combination is also valid.

## Validity

Identical to v14: every admissible $f$ produces a feasible point at
$\Omega^{*} = \|f * f\|_\infty$.  Hence
$\Omega^{(16)} \le C_{6.2}$.

## Numerical observation

At $N = 8$, $K = 8$, MOSEK returns $\Omega^{(16)} = 1.033$ —
**identical** to v14 at the same parameters.  Adding the Bochner
Toeplitz and (W) constraints does not improve the bound at this scale.

The reason is that the binding constraint is the autocorrelation
Fejér, and the Bochner Toeplitz / (W) constraints are slack at the
LP optimum.  This suggests that the bottleneck is the linear bounds
on $a_k, b_k$ from $p$ (which the LP can satisfy by choosing a
non-uniform $p$ that minimizes $\sum |\hat f(k)|^{4}$), not the
PSD/(W)/Bochner constraints.

To improve the bound, we need either:
1. Tighter $a_k, b_k$ bounds (e.g., higher-order f-side moment
   matching).
2. A constraint that **forces** $\sum |\hat f(k)|^{2}$ from below
   (Parseval), which is hard at finite $K$ due to tail control.
3. Larger $N$ and $K$ (slow but reliable).

v17 will explore option 3 (larger parameters) with v14 alone to see
how high the bound can go.
