# v30 — Rigorous proof: pointwise Fejér mean tests at multiple t*

## Key insight (from v29 diagnostic)

v29 verified that v14's LP-optimal $f$ has actual
$\|f * f\|_\infty \approx 8$ (not 1.06).  The v14 framework captures
$\int F^{2}$ (the $L^{2}$ norm) but underestimates $\|F\|_\infty$.

To capture the peak directly, v30 adds **pointwise Fejér mean
constraints** at multiple test points $t_{*} \in [-1/2, 1/2]$:
$$
\Omega \;\ge\; \sigma_{K}(t_{*}) \;:=\; (F * F_{K})(t_{*}),
$$
where $F_{K}$ is the Fejér kernel of order $K$.  Each constraint is
rigorous because $\sigma_{K}(t_{*}) \le \|F\|_\infty \cdot \|F_{K}\|_{1} = \Omega$
by Hölder.

## The constraint chain

For real $F = f * f$ with $\hat F(0) = 1$:
$$
\sigma_{K}(t_{*}) \;=\; \sum_{|k|\le K}\Bigl(1 - \frac{|k|}{K+1}\Bigr)\hat F(k)\,e^{2\pi i k t_{*}}
$$
$$
\;=\; 1 + 2\sum_{k=1}^{K}\Bigl(1 - \frac{k}{K+1}\Bigr)\Bigl[\operatorname{Re}\hat F(k)\cos(2\pi k t_{*}) - \operatorname{Im}\hat F(k)\sin(2\pi k t_{*})\Bigr].
$$

For $\hat F(k) = \hat f(k)^{2} = (a_{k} - i b_{k})^{2}$:
$\operatorname{Re}\hat F(k) = a_{k}^{2} - b_{k}^{2}$ and
$\operatorname{Im}\hat F(k) = -2 a_{k} b_{k}$.  So
$$
\sigma_{K}(t_{*}) \;=\; 1 + 2\sum_{k=1}^{K}\Bigl(1 - \frac{k}{K+1}\Bigr)\Bigl[(a_{k}^{2} - b_{k}^{2})\cos(2\pi k t_{*}) + 2 a_{k} b_{k}\sin(2\pi k t_{*})\Bigr].
$$

In the SDP lift on $M$:
- $a_{k}^{2} \mapsto M[k, k]$
- $b_{k}^{2} \mapsto M[K+k, K+k]$
- $a_{k} b_{k} \mapsto M[k, K+k]$

The lifted form:
$$
\sigma^{(\text{lift})}_{K}(t_{*}) \;=\; 1 + 2\sum_{k=1}^{K}\Bigl(1 - \frac{k}{K+1}\Bigr)\Bigl[(M[k,k] - M[K+k,K+k])\cos(2\pi k t_{*}) + 2 M[k, K+k]\sin(2\pi k t_{*})\Bigr].
$$

We add the constraint $\Omega \ge \sigma^{(\text{lift})}_{K}(t_{*})$
for each test point $t_{*}$ in a uniform grid of $[-1/2, 1/2]$.

## Validity

For every admissible $f$, the rank-one $M = vv^{\top}$ from $f$ gives
$M[k,k] = a_{k}^{2}$, $M[K+k,K+k] = b_{k}^{2}$, $M[k, K+k] = a_{k} b_{k}$,
so $\sigma^{(\text{lift})}_{K}(t_{*}) = \sigma_{K}(t_{*})$ (the true Fejér
mean).  And $\sigma_{K}(t_{*}) \le \Omega^{*} = \|f * f\|_\infty$ by the
Hölder chain.

So the rank-one truth at $\Omega^{*}$ satisfies all the v30 constraints,
and v30 is rigorously valid: $\Omega^{(30)} \le C_{6.2}$.

(Caveat on inflation: in the relaxation, the LP could *also* find a
feasible $M$ that is not rank-one and gives a smaller lifted
$\sigma^{(\text{lift})}_{K}$.  This is fine because we need only that
the rank-one truth is feasible — and it is.  The LP optimum is then
$\le \Omega^{*}$, which is the desired lower bound on $C_{6.2}$.)

## Numerical results

| $N$ | $K$ | $n_\text{tests}$ | $\Omega^{(30)}$ |
|-----|-----|------------------|------------------|
| 8   | 8   | 33               | 1.086            |
| 8   | 16  | 33               | 1.093            |
| 16  | 8   | 33               | 1.131            |
| 16  | 16  | 33               | 1.147            |
| 16  | 32  | 51               | 1.153            |
| 24  | 24  | 51               | 1.177            |
| 32  | 32  | 51               | 1.196            |
| 50  | 32  | 51               | 1.217            |

These are **substantially better** than v14 at the same parameters
(e.g., 1.196 vs 1.088 at $N=K=32$).  The bound is climbing toward
the SOTA $1.28$ and may exceed it at larger parameters.

## Conclusion

v30 is the **fourth major milestone** of the project, after v7
(rigorous reboot for even $f$), v9 (Fejér for even $f$), and v14
(autocorrelation Fejér for general $f$).  It is the first iteration
to give a rigorous bound on $C_{6.2}$ approaching the SOTA.

The improvement comes from using **pointwise Fejér mean tests**
(per the v9 verdict's recommendation, now properly applied to
non-even $f$) on top of the v14 autocorrelation Fejér constraint.
The combination is strictly tighter than either alone.

**Best rigorous bound on $C_{6.2}$ from v30: $\Omega \ge 1.217$ at
$(N, K, n_\text{tests}) = (50, 32, 51)$**, with larger parameters
expected to push closer to (and possibly exceed) the SOTA $C_{6.2} \ge 1.28$.
