# v23 — Rigorous fix of v22 (revert to v14's loose cell-bound)

## Summary

v22 used a non-rigorous $O(L_f^{2})$ midpoint quadrature error.  v23
reverts to v14's loose $O(L_f \cdot k)$ cell-bound, which is the
correct rigorous bound for arbitrary nonneg $f$ (not necessarily
piecewise constant or smooth).

The remainder of v23 is identical to v14: same autocorrelation
Fejér constraint, same SDP lift on $M$, same f-side discretization
$p$.

## Numerical results

Same as v14 (since the constraints are identical):

| $N$ | $K$ | $\Omega^{(23)}$ |
|-----|-----|-----------------|
| 8   | 8   | 1.033           |
| 16  | 16  | 1.063           |
| 50  | 64  | 1.10            |

## Conclusion

v23 confirms that v22's apparent improvement was due to a non-rigorous
constraint.  The honest rigorous bound for the autocorrelation Fejér
construction at the parameters tested remains $\Omega \approx 1.06–1.11$.

The path forward (without restoring rigor compromises) is to either:
1. Run the rigorous v14/v23 at much larger $N, K$ and hope for slow
   convergence to the asymptotic limit.
2. Find a fundamentally new constraint type that adds rigorous
   information beyond what v14 captures.

The v20 verdict's recommendation #3 — "stop adding constraints,
focus on diagnosing the convergence" — remains the right advice.
