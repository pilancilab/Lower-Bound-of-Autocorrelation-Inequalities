# Verdict — Autocorrelation Inequality 6.2, v28

## 1. Summary

v28 adds the **Bochner Toeplitz PSD constraint on \(\hat F\)** —
exactly the dual-side Bochner constraint I recommended in the v13
verdict's "v14 plan" (Comment 4).

The constraint: \(T_F = [\hat F(j-i)]_{i,j=0}^{K} \succeq 0\) (Hermitian
PSD), with \(\hat F(k) = \hat f(k)^{2}\) lifted via the moment matrix
\(M\) as
\[
\operatorname{Re}\hat F(k) = M[k,k] - M[K+k, K+k],\quad \operatorname{Im}\hat F(k) = -2 M[k, K+k].
\]

The constraint is mathematically valid (Bochner's theorem for the
positive measure \(F = f*f\)), but produces *no numerical
improvement* over v14: the LP optimum is identical at every tested
\((N, K)\).

The proposer correctly identifies this as the *last* constraint
addition and explicitly declares the v14 ceiling.

## 2. Verdict

**VALID** (rigorous, no improvement, but the right analysis).

### Verification

**Bochner's theorem**: For \(F = f*f\) with \(f\ge 0\), the Fourier
sequence \((\hat F(k))_{k\in\mathbb Z}\) is positive definite (i.e.,
the Hermitian Toeplitz \([\hat F(j-i)]\) is PSD at every size). ✓

**Validity in the lift**: at the rank-one lift \(M^{\star} = vv^{\top}\),
the lifted entries match the actual \(\hat F(k) = \hat f(k)^{2}\),
and Bochner gives \(T_F\succeq 0\). For non-rank-one \(M\) in the
relaxation, adding \(T_F\succeq 0\) restricts the feasible set
further, preserving validity (every rank-one lift is still feasible).

### Why it doesn't help

The proposer's observation:
> v14's LP-optimal \(M\) already has lifted \(T_F\) PSD (because the
> LP uses the rank-one truth \(M = vv^{\top}\) via the
> \(u_k = M[k,k]+M[K+k,K+k]\) construction, which forces \(T_F\) to
> be PSD automatically).

is correct *at the LP optimum*. v14's SOC chain
\(v_k \ge u_k^{2}\) effectively forces the LP to behave as if
\(M\) were rank-one at the optimum, and the rank-one lift always
satisfies Bochner's constraint. So adding the constraint doesn't
exclude the LP's optimal point.

In a more relaxed (non-rank-one) \(M\), the constraint might
exclude some points — but those points don't appear at the LP
optimum.

## 3. The "v14 ceiling" declaration

The proposer explicitly declares:
> This is the **last constraint type** I attempted to add to v14
> without sacrificing rigor. All Fourier-side and PSD-side
> restrictions have been tried; none improves the bound.

After 28 iterations on 6.2, the project's best rigorous bounds are:

| Constant | Best bound | Iteration |
|---|---|---|
| \(C_{6.2}^{\text{even}}\) | 1.689 | v9 (\(N=50, K=64\)) |
| \(C_{6.2}\) (general) | **1.111** | v17/v25 (\(N=100, K=32\)) |
| Asymptotic limit of v14 | \(\sim 1.07\)–\(1.15\) | v26 diagnostic |
| SOTA on \(C_{6.2}\) | 1.28 | (unchanged) |

**The v14 framework cannot reach the SOTA**, and after 14 iterations
of attempted tightenings (v15–v28), no new constraint type produces
any improvement.

The proposer's correct conclusion: the path forward requires either
(a) running v14 at much larger parameters (probably futile given the
\(\sim 1.07\)–\(1.15\) asymptotic ceiling identified in v26), or
(b) a fundamentally different relaxation framework (Lasserre on
moments of \(f\), SOS hierarchy, dual constructions, etc.).

## 4. Comments

### Comment 1. The 6.2 project has effectively concluded with this iteration.

After v28, the proposer has tried all the standard constraint types
in the v14 framework:

| Iteration | Constraint type | Result |
|---|---|---|
| v15 | shifted (W) windows | no improvement |
| v16 | v13 Bochner on \(\hat f\) | no improvement |
| v18 | L³ constraint | no improvement |
| v19 | multi-\(K\) Fejér | tiny improvement |
| v20 | Lasserre level-2 on \((a, b)\) | no improvement |
| v21 | squared Fejér kernel | tiny improvement |
| v24 | \(M_\infty\) variable | no improvement |
| **v28** | **Bochner Toeplitz on \(\hat F\)** | **no improvement** |

Plus the diagnostic iterations (v26, v27) that identified the
fundamental ceiling.

This is a *thorough* exploration of the constraint space within the
v14 framework, and it reaches the right negative conclusion: v14
cannot reach the SOTA via *any* of these tightenings.

### Comment 2. The v28 / v13 dual.

In the 6.2 v13 verdict's Comment 1 ("v13 is necessary but not
sufficient"), I noted that the f-side Bochner constraint (v13) only
constrains \((a_k, b_k)\), not the *higher-order* moments
\(M[k,k]\), \(M[K+k, K+k]\), \(M[k, K+k]\) in the lifted SDP.

v28 adds the *F-side* Bochner constraint on the lifted moments of
\(\hat F\), which directly involves the diagonal entries
\(M[k,k]\), \(M[K+k, K+k]\). This *should* in principle constrain
the diagonal inflation that defeated v10 / v12 / v13. 

However, by the time v28 is added on top of v14, the v14 SOC chain
already eliminates the diagonal inflation (because the LP optimum
is at the rank-one boundary anyway). So v28's constraint, while
*technically* what was needed in v13, is redundant on top of v14.

This is a subtle structural observation: **the right constraint
needs to be added at the right place in the relaxation hierarchy.**
Adding v28 to v14 is too late to break new ground.

### Comment 3. Inheritance is clean.
v28 inherits from v14 (no v0 hand-wave), and the new constraint is
internally consistent.

## 5. Final feedback to the proposer

After 28 iterations, the 6.2 project has:

1. **Established a rigorous lower bound** \(C_{6.2}\ge 1.111\) on
   the unrestricted constant via v14/v17's autocorrelation Fejér
   construction.
2. **Diagnosed the v14 ceiling** at \(\sim 1.07\)–\(1.15\) (v26),
   below the SOTA 1.28.
3. **Exhausted the standard constraint additions** to v14 (v15,
   v16, v18, v19, v20, v21, v24, v28) without improvement.

The project has produced *real* mathematical content (the v14
construction is novel and rigorous), but has not reached the SOTA.
The honest path forward is to either:

* **Declare v25 as the final result** (\(C_{6.2}\ge 1.111\)) and
  document the project's contribution as a "first non-even rigorous
  bound on \(C_{6.2}\) via convex relaxation," noting that it's
  below the SOTA.
* **Or pivot to a fundamentally different framework** for v29+,
  acknowledging that the v14 family has been exhausted.

I would recommend the former. The v25 synthesis is a clean
declaration, and pushing further on v14 is unlikely to help.

---

**Final verdict: VALID** (rigorous, no improvement; the v14 ceiling
declaration is honest and correct).
