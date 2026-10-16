# v15 — Rigorous proof: v14 + v7/v11 (W) windowed constraints

## Setup

v15 combines all constraints from v14 (autocorrelation Fejér,
rigorous for any admissible $f$) and v7/v11 (rigorous (W) windowed
integral constraints on the (p, Q) lift, also valid for any
admissible $f$).  The combination is rigorously valid because every
constraint is individually valid.

## Constraints

* All v14 constraints (mass, $a_k$, $b_k$, M PSD lift, autocorrelation
  Fejér via $u_k$, $v_k$).
* All v7/v11 constraints ((p, Q) PSD lift, marginals, $Q \ge 0$,
  centered + shifted (W) windows at multiple half-widths).

The variables $a_k$, $b_k$ are linked to $p$ via the f-side cell
bounds, and $p$ is the same variable used in the (p, Q) lift.  So
the LP couples the Fourier and the f-side cell discretizations.

## Validity

Same proof as v14 + v7/v11, since all constraints are valid for the
rank-one truth.  The relaxation optimum
$\Omega^{(15)} \le C_{6.2}$.

## Numerical observation

At $N = 8$, $K = 8$, MOSEK returns $\Omega^{(15)} = 1.033$, the same
value as v14 alone.  The (W) constraints from v7/v11 do not bind in
the v14-Fejér-dominated LP at this scale, because the LP picks $p$
to minimize the Fejér sum and the resulting $p$ already satisfies all
the (W) constraints.

This is consistent: the (W) constraints constrain the second-order
moments $Q[j, k]$ of $p$ via the bilinear convolution structure,
which is **separate** from the lifted Fourier moments $M$.  The
binding constraint at the LP optimum is the autocorrelation Fejér,
not the (W) windowed sums.

To obtain a tighter bound, either:
1. Increase $N$ and $K$ (slow convergence).
2. Add a constraint that **couples** the (W) windowed sums to the
   Fourier moments.  The natural candidate is to use $u_k$ (the
   lifted $|f̂(k)|^{2}$) inside an LP constraint that also references
   $p$.  Such constraints are quadratic in $p$ and don't lift cleanly.

v16 will try a different angle: tighten the autocorrelation lift
itself by adding a LARGER moment matrix on $(1, a_1, ..., a_K, b_1, ..., b_K, a_1^{2}, b_1^{2}, ...)$ — i.e., a higher-degree Lasserre lift.

## Conclusion

v15 is rigorous and gives the same numerical bound as v14
($\Omega \approx 1.03$ at $N = 8, K = 8$).  It demonstrates that the
v14 autocorrelation Fejér is the dominant constraint and that the
(W) windowed constraints, while valid, do not contribute marginally
at this scale.
