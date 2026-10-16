# Verdict — Autocorrelation Inequality 6.2, v21

## 1. Summary

v21 follows the v19 verdict's recommendation to try a different
positive trigonometric polynomial in the v14 framework. The chosen
kernel is the **squared Fejér kernel**
\[
P(t) := \frac{F_K(t)^{2}}{Z_K},\qquad Z_K := \int F_K^{2}\,dt.
\]
The constraint is the v14 Hölder chain with \(P\) replacing \(F_K\):
\(\Omega\ge 1 + 2\sum_{m=1}^{2K}\hat P(m)|\hat F(m)|^{2}\).

The numerical bound at \((N, K) = (16, 16)\) is \(\Omega^{(21)} = 1.068\),
compared to v14's 1.063 at the same parameters. A marginal
improvement of \(\sim 0.005\).

## 2. Verdict

**VALID.** The squared Fejér kernel is a valid positive
trigonometric polynomial, and the v14 Hölder chain works for any
such kernel. The math is correct.

## 3. Verification

### Property (P1): \(P\ge 0\). ✓
\(P(t) = F_K(t)^{2}/Z_K\) is the square of a real number divided by
a positive constant. Nonneg. ✓

### Property (P2): \(\int P = 1\). ✓
By normalization \(Z_K = \int F_K^{2}\,dt\). ✓

### Property (P3): \(\hat P(m) \ge 0\). ✓
\(P = F_K^{2}/Z_K\) has Fourier coefficients
\(\hat P(m) = \widehat{F_K^{2}}(m)/Z_K\). The Fourier transform of
the product \(F_K \cdot F_K\) is the convolution
\(\hat F_K * \hat F_K\):
\[
\widehat{F_K^{2}}(m) = (\hat F_K * \hat F_K)(m) = \sum_{k}\hat F_K(k)\hat F_K(m-k).
\]
Since \(\hat F_K(k) = (1 - |k|/(K+1))_{+}\ge 0\) for all \(k\), the
convolution is a sum of nonneg terms, hence nonneg. So
\(\hat P(m)\ge 0\). ✓

(This is also visible from \(P\ge 0\) plus Bochner's theorem in
reverse: nonneg \(P\) implies positive-definite Fourier sequence. The
convolution computation gives the same result more directly.)

### The constraint chain. ✓

Same as v14 with \(P\) replacing \(F_K\):
\[
\int P\,(F * \check F)\,dt \;\le\; \|P\|_{1}\cdot\|F*\check F\|_{\infty} \;=\; 1\cdot \int F^{2} \;\le\; \Omega\cdot \int F \;=\; \Omega.
\]
* \(\|P\|_{1} = \int P = 1\) by (P2).
* \(\|F*\check F\|_{\infty} = (F*\check F)(0) = \int F^{2}\) (v14, Property C).
* \(\int F^{2}\le \Omega\int F = \Omega\) (using \(F\le \Omega\) and \(\int F = 1\)).

✓ All steps correct.

### Parseval expansion. ✓
\[
\int P(F*\check F)\,dt = \sum_{m\in\mathbb Z}\hat P(m)\overline{\widehat{F*\check F}(m)} = \sum_{m}\hat P(m)|\hat F(m)|^{2},
\]
using \(\widehat{F*\check F}(m) = |\hat F(m)|^{2}\) (real and nonneg).

The sum is supported on \(|m|\le 2K\) because \(P = F_K^{2}\) is a
trigonometric polynomial of degree \(2K\) (since \(F_K\) has degree
\(K\) and the product has degree \(2K\)). So:
\[
\sum_{m\in\mathbb Z}\hat P(m)|\hat F(m)|^{2} = \hat P(0) + 2\sum_{m=1}^{2K}\hat P(m)|\hat F(m)|^{2} = 1 + 2\sum_{m=1}^{2K}\hat P(m)|\hat F(m)|^{2}.
\]
(Using \(\hat P(0) = 1\) since \(\int P = 1\) and the symmetry
\(|\hat F(-m)| = |\hat F(m)|\) for real \(F\).) ✓

### The lifted constraint. ✓
Using the v14 lift \(v_m\ge u_m^{2}\) with
\(u_m = M[m,m] + M[K+m, K+m]\) and \(M[m,m]\ge a_m^{2}\),
\(M[K+m, K+m]\ge b_m^{2}\):
\[
v_m \;\ge\; u_m^{2} \;\ge\; (a_m^{2}+b_m^{2})^{2} \;=\; |\hat f(m)|^{4} \;=\; |\hat F(m)|^{2}.
\]
So replacing \(|\hat F(m)|^{2}\) by \(v_m\) in the constraint
gives a *stronger* constraint:
\[
\Omega \ge 1 + 2\sum_{m=1}^{2K}\hat P(m)\,v_m.
\]
And the rank-one lift achieves equality, so the constraint is valid
for the rank-one lift of any admissible \(f\). ✓

(Note: the LP needs \((a_m, b_m)\) for \(m\) up to \(2K\), not just
\(K\). The "effective degree" is \(2K\), as the proposer notes in
the numerical table.)

## 4. Numerical comparison

| Iteration | \((N, K)\) | \(\Omega\) | Effective degree |
|---|---|---|---|
| v14 | \((16, 16)\) | 1.063 | 16 |
| v14 | \((16, 32)\) | 1.066 | 32 |
| **v21** | **\((16, 16)\)** | **1.068** | **32** |

So at the same *effective degree* (32), v21 beats v14 by \(\sim 0.002\).
A marginal improvement, but it shows that *kernel choice matters*:
the squared Fejér is a slightly better test polynomial than the
plain Fejér at the same degree.

## 5. Comments

### Comment 1. The improvement is in the right direction.

Adding more positive-trig-polynomial test functions to the v14
framework is the right move (per my v19 verdict §5 Comment 1).
v21 demonstrates that different kernels give *different* (slightly
better or worse) bounds, and the optimal choice is non-trivial.

### Comment 2. The optimal positive trig polynomial.

For *fixed degree* \(K\), the *optimal* nonneg trigonometric
polynomial \(P\) of degree \(\le K\) — in the sense of giving the
*largest* LP bound when used as a test function — is itself the
solution to a sub-optimization problem. The optimal \(P\) depends
on the LP variables (which is what makes it non-trivial) and is in
general computed via the *Carathéodory–Toeplitz parameterization*:
\(P\) is parameterized by an auxiliary PSD Toeplitz matrix.

I would suggest in v22 to *jointly optimize* the test function \(P\)
and the LP variables \((p, a, b, M)\) by adding the
Carathéodory–Toeplitz parameterization. This is a standard SDP
reformulation, and it gives the *best possible* bound from the
positive-trig-polynomial test approach at any fixed degree.

### Comment 3. Convergence is still slow.

v21 gives 1.068 at \((N, K) = (16, 16)\). v17 gives 1.111 at
\((N, K) = (100, 32)\). The plateau around 1.10–1.15 is unchanged.
The squared Fejér kernel does not break it.

For the project to actually beat the SOTA \(C_{6.2}\ge 1.28\), the
proposer needs *qualitatively new* ideas, not just minor variations
on the v14 theme. As I noted in v19 §5 Comment 3, the binding
bottleneck is either (a) the discretization on \(p\to a_k\), or
(b) the actual infimum of \(\int F^{2}\) being smaller than \(4/3\).
The proposer should diagnose which.

### Comment 4. Inheritance is clean.
v21 inherits from v14 (no v0 hand-wave). The construction is clean.

## 6. Feedback to the proposer

The squared Fejér kernel is a good first move on the
"different test polynomial" axis. A few more directions to try:

1. **Jackson kernel**: \(J_K(t) = c_K F_{\lfloor K/2\rfloor}(t)^{2}\)
   normalized to \(\int J_K = 1\). Same idea as squared Fejér but
   with different scaling, gives different rates.
2. **de la Vallée–Poussin kernel**: \(V_K = 2 F_{2K-1} - F_{K-1}\).
   This is a *positive* trig polynomial that is constant 1 in
   Fourier space for \(|k|\le K\), then linearly tapers. Different
   weight profile.
3. **Optimal positive trig polynomial via Carathéodory–Toeplitz.**
   Jointly optimize the test function and the LP variables. This is
   the strongest version and would be the *best* result obtainable
   from positive-trig-polynomial methods at any fixed degree.

If v22 implements (3), I expect a real numerical improvement over
v21 — possibly enough to break the 1.10 plateau.

---

**Final verdict: VALID** (marginal improvement; right direction).
