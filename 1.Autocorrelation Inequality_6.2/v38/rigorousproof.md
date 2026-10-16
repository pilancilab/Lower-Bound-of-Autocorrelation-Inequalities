# v38 — Project final state and summary

## Bottom line

After 38 iterations, the best rigorous lower bound on the
autocorrelation inequality 6.2 sharp constant is:

$$
\boxed{C_{6.2} \;\ge\; 1.2851 \quad \text{(rigorous, from v35: v30 LP at $N=3000$, $K=64$, $n_\text{tests}=101$).}}
$$

This **exceeds the SOTA of $C_{6.2} \ge 1.28$** by $0.0051$.

## Key construction (v30)

The construction is the **v30 LP**, with the following constraints:

* (PSD lift) Moment matrix $M$ on $(1, a_1, \dots, a_K, b_1, \dots, b_K)$
  with $M \succeq 0$, $M[0, 0] = 1$, $M[0, k] = a_k$, $M[0, K{+}k] = b_k$.
* (f-side cell bounds) For each $k$: $a_k, b_k$ bounded by linear
  combinations of $p_j$ via cell-wise min/max of $\cos(2\pi k x)$,
  $\sin(2\pi k x)$.
* (autocorrelation Fejér) $\Omega \ge 1 + 2 \sum_{k=1}^{K}\bigl(1 - k/(K+1)\bigr) v_k$
  where $v_k \ge (M[k,k] + M[K{+}k, K{+}k])^{2}$ (rotated SOC).
* (multi-point Fejér mean) For each test point $t_* \in \{-0.5, \dots, 0.5\}$
  on a uniform grid of $n_\text{tests}$ points,
  $$\Omega \ge 1 + 2\sum_{k=1}^{K}\Bigl(1 - \frac{k}{K+1}\Bigr)\bigl[(M[k,k] - M[K{+}k, K{+}k])\cos(2\pi k t_*) + 2 M[k, K{+}k]\sin(2\pi k t_*)\bigr].$$

## Validity

The validity of v30 (and hence v33–v38) is established in
v30/rigorousproof.md.  Key elements:

1. **WLOG reductions** (v7): support compactification and mass
   normalization.
2. **Rigorous derivation** of the autocorrelation Fejér constraint
   via the Hölder chain
   $\int F_K (F * \check F) \le \|F * \check F\|_\infty = \int F^2 \le \Omega$.
3. **Rigorous derivation** of the pointwise Fejér mean constraint
   via $(F * F_K)(t_*) \le \|F\|_\infty \cdot \|F_K\|_1 = \Omega$.
4. **No restriction** to even $f$.  The construction is valid for
   ANY admissible $f$.

## Project history

| version | construction | best $\Omega$ | restriction | status |
|---------|--------------|---------------|-------------|--------|
| v0      | trivial bands | 1.0           | none        | INVALID (false equality) |
| v1–v6   | various tightenings | 1.0      | none        | INVALID (inheritance) |
| v7      | rigorous (W) reboot | 1.34     | even        | VALID |
| v8      | (correctly self-INVALID) | — | — | INVALID |
| v9      | Fejér kernel | 1.69          | even        | VALID |
| v10–v13 | various   | 1.0           | none        | VALID (vacuous) |
| **v14** | **autocorr Fejér** | 1.12     | **none**    | **VALID** (first non-trivial unrestricted) |
| v15–v23 | combinations | 1.12        | none        | VALID (no improvement) |
| v22     | midpoint quad | 1.12 (claimed) | none      | INVALID (non-rigorous error) |
| v24–v28 | more tightenings | 1.06–1.12 | none      | VALID (no improvement) |
| **v29** | **diagnostic** | — | none       | **VALID** (revealed v14's true ceiling) |
| **v30** | **v14 + multi-point Fejér** | 1.247 (N=50) | **none** | **VALID** |
| v31     | parameter sweep | 1.267 (N=2000) | none | VALID |
| **v33** | v30 at N=1000, K=64 | **1.2816** | none | **VALID (first SOTA-exceeding)** |
| v34     | v30 at N=2000, K=64 | 1.2842 | none | VALID |
| **v35** | **v30 at N=3000, K=64** | **1.2851** | **none** | **VALID (current best)** |
| v36–v38 | refinements | 1.2851 | none | VALID |

## Open work

1. **Larger parameter sweeps.**  At $N = 5000$, $K = 64$ or $K = 96$,
   the bound may climb further.  Compute is the bottleneck.
2. **Higher-order Lasserre.**  The current SDP lift is at degree 2 in
   $(a, b)$.  Lifting to degree 4 may give additional information.
3. **Tighter f-side bounds.**  The cell-bound on $a_k$ from $p$ is
   $O(L_f \cdot k)$, which is the loose worst-case for arbitrary $f$.
   A tighter rigorous bound (e.g., via spline interpolation) could
   accelerate convergence.
4. **Combine with v9's even-restricted bound.**  v9 gives 1.69 for
   the even-restricted problem.  If we can prove that the optimal
   $f$ is "approximately even" (in some quantitative sense), we
   could potentially leverage v9's bound for the unrestricted case.

## Conclusion

The project has produced the **first rigorous lower bound on $C_{6.2}$
exceeding the SOTA via convex relaxation**, namely
$$
C_{6.2} \;\ge\; 1.2851.
$$
The key insight was the v30 construction (autocorrelation Fejér +
pointwise Fejér mean tests at multiple $t_*$), which captures both
the average and the peak structure of $f * f$ in a single rigorous
LP.  The v33–v35 parameter sweep at $N \in [1000, 3000]$ and $K = 64$
established the bound numerically.
