# v20 — Rigorous proof: Lasserre level-2 lift on $(a, b)$

## Setup

v20 strengthens the v14 SDP lift on $M$ from degree 1 to degree 2 in
the variables $(a_1, \dots, a_K, b_1, \dots, b_K)$.  The level-2
moment matrix has rows/cols indexed by all monomials of degree
$\le 2$, and is constrained PSD.

## Constraints

* All v14 constraints on $(p, a, b, M)$ (Fourier f-side bounds,
  PSD lift on $M$, autocorrelation Fejér via $v_k \ge u_k^{2}$).
* New: a level-2 moment matrix $M_{2}$ on the $2K$ variables
  $(a, b)$, indexed by multi-indices of degree $\le 2$, constrained
  $M_{2} \succeq 0$.

The identifications:
* $y_{0} = 1$
* $y_{e_{k}^{a}} = a_{k}$, $y_{e_{k}^{b}} = b_{k}$
* $y_{2e_{k}^{a}} = a_{k}^{2}$, $y_{2e_{k}^{b}} = b_{k}^{2}$,
  $y_{e_{k}^{a} + e_{k}^{b}} = a_{k} b_{k}$, etc.

## Validity

Every admissible $f$ produces a feasible point of v20: the rank-one
moment sequence $y_{\alpha} = a^{\alpha_a} b^{\alpha_b}$ has the
moment matrix $M_{2}$ rank-one PSD.  Hence v20 is rigorously valid.

## Numerical observation

| $N$ | $K$ | $\Omega^{(20)}$ |
|-----|-----|-----------------|
| 8   | 3   | 1.022           |
| 8   | 4   | 1.026           |
| 16  | 3   | 1.028           |
| 16  | 4   | 1.040           |

These match v14 at the same parameters.  The level-2 lift does not
tighten the bound at this scale.

The reason is that the LP optimum at v14 already corresponds to a
rank-one moment matrix on $(1, a, b)$, and lifting to higher-degree
moments simply adds more rank-one PSD constraints that are
automatically satisfied.

To get strictly more, we would need to add **additional algebraic
constraints** that bind only with the higher-degree moments — for
example, the explicit Parseval identity at finite $K$, or specific
sum-rules from the autocorrelation Toeplitz structure on $\hat F$.

## Conclusion

v20 is rigorous but produces no improvement over v14.  Higher-order
SDP lifts on $(a, b)$ alone do not help; the binding constraint
remains the v14 autocorrelation Fejér plus the linear bounds from
$p$.
