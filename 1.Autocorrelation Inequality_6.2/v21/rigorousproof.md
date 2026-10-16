# v21 — Rigorous proof: squared Fejér kernel test

## Setup

v21 follows the v19 verdict's recommendation to try a different
positive trigonometric polynomial.  We use the **squared Fejér
kernel**:
$$
P(t) := \frac{F_K(t)^{2}}{Z_K},
\qquad Z_K := \int F_K(t)^{2}\,dt = \sum_{|k|\le K}\Bigl(1 - \frac{|k|}{K+1}\Bigr)^{2}.
$$

**Key properties.**

(P1) $P(t) \ge 0$ pointwise (it is the square of a real function).

(P2) $\int P\,dt = 1$ by normalization.

(P3) $P$ is a trigonometric polynomial of degree $2K$, so its
Fourier coefficients are supported on $|m| \le 2K$.  Specifically,
$$
\hat P(m) \;=\; \frac{1}{Z_K}\sum_{k} F_K^{\,\hat{}\,}(k)\,F_K^{\,\hat{}\,}(k+m)
\;=\; \frac{1}{Z_K}\sum_{k}\Bigl(1 - \frac{|k|}{K+1}\Bigr)\Bigl(1 - \frac{|k+m|}{K+1}\Bigr).
$$
This is the **discrete autocorrelation** of the Fejér Fourier
coefficients.  All $\hat P(m) \ge 0$ (since $P \ge 0$).

## The constraint

By the same Hölder + L²-bound chain as v14:
$$
\int P(t)\,(F * \check F)(t)\,dt \;\le\; \|P\|_1 \cdot \|F * \check F\|_\infty \;=\; 1 \cdot \int F^{2} \;\le\; \Omega.
$$
The left-hand side, by Plancherel,
$$
\int P\,(F*\check F)\,dt
\;=\; \sum_{|m|\le 2K} \hat P(m)\,|\hat F(m)|^{2}
\;=\; 1 + 2\sum_{m=1}^{2K} \hat P(m)\,|\hat F(m)|^{2},
$$
using $\hat P(0) = 1$ and the symmetry $|\hat F(-m)| = |\hat F(m)|$.

So
$$
\Omega \;\ge\; 1 + 2\sum_{m=1}^{2K} \hat P(m)\,|\hat F(m)|^{2}.
\tag{$\star_{21}$}
$$

The lifted form (using the same $v_m \ge u_m^{2}$ chain as v14):
$$
\Omega \;\ge\; 1 + 2\sum_{m=1}^{2K} \hat P(m)\,v_m,
$$
which is a single linear constraint on $(\Omega, v_1, \dots, v_{2K})$.

## Validity

Same as v14: the rank-one truth $M = vv^{\top}$ from any admissible
$f$ gives $u_m = a_m^2 + b_m^2 = |\hat f(m)|^{2}$ and
$v_m = u_m^{2} = |\hat F(m)|^{2}$, so the constraint matches
$(\star_{21})$ at equality.  Hence $\Omega^{(21)} \le C_{6.2}$.

## Numerical results

| $N$ | $K$ (eff degree $2K$) | $\Omega^{(21)}$ |
|-----|------------------------|-----------------|
| 8   | 4 (8)                  | 1.034           |
| 8   | 8 (16)                 | 1.038           |
| 8   | 16 (32)                | 1.039           |
| 16  | 4 (8)                  | 1.057           |
| 16  | 8 (16)                 | 1.065           |
| 16  | 16 (32)                | 1.068           |

For comparison, v14 at $(N, K) = (16, 16)$ gives $\Omega = 1.063$,
and v21 at the same parameters gives $\Omega = 1.068$ (a 0.005
improvement).  At $(16, K=16)$, the squared Fejér uses an effective
degree of 32, so the comparison is more fair to v14 at $K = 32$
which gave $\Omega = 1.066$ (so v21 still beats by 0.002).

The improvement is real but small.  Squared Fejér is *one* family
of positive trig polynomials; the optimal one (over all positive
trig polynomials of bounded degree) would give the largest LP
bound, but is itself a non-trivial optimization.

## Conclusion

v21 demonstrates that different positive trig polynomial test
functions give marginally different bounds in the v14 framework.
The squared Fejér is slightly tighter than the regular Fejér for
the same number of frequencies.  Future work could explore optimal
test polynomials.
