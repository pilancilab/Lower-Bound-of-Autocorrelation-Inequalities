# v13 — Rigorous proof: v12 + Bochner Toeplitz, observed failure

## Setup

v13 inherits all v12 constraints and adds the **Bochner / Toeplitz**
PSD constraint on the (real, imaginary) parts of $\hat f(k)$:

$$
T \;:=\; \bigl[\hat f(j-k)\bigr]_{j,k=0}^{K}
\;\succeq\; 0\quad\text{(Hermitian PSD)},
$$

where $\hat f(0) = 1$, $\hat f(k) = a_k - i b_k$ for $k > 0$, and
$\hat f(-k) = \overline{\hat f(k)} = a_k + i b_k$.  In real form:
$$
T \;=\; T_R + i T_I,\qquad
T_R[i,j] = a_{|i-j|},\;\;
T_I[i,j] = \mathrm{sign}(j-i)\,b_{|i-j|}.
$$
The Hermitian PSD condition is equivalent to the real-block PSD
condition
$$
\begin{pmatrix} T_R & -T_I \\ T_I & T_R \end{pmatrix} \;\succeq\; 0.
$$

## Why this is rigorous

**Bochner's theorem.**  A complex sequence $(c_k)_{k\in\mathbb Z}$ is
positive definite (i.e., the Hermitian Toeplitz matrices
$[c_{i-j}]_{i,j=0..N}$ are PSD for every $N$) if and only if there
exists a nonneg Borel measure $\mu$ on $[-1/2, 1/2]$ such that
$c_k = \int e^{-2\pi i k x}\,d\mu(x)$.

For our admissible $f \ge 0$, the Fourier coefficient sequence
$\hat f(k)$ is by definition a positive-definite sequence, and the
Toeplitz matrix $T$ is PSD.  Hence the constraint $T \succeq 0$
is **valid** for every admissible $f$, and adding it to the
relaxation only restricts the feasible set in a consistent way.

## Why v13 still fails numerically

At $N \in \{8, 12\}$ and $K \in \{6, 12\}$, MOSEK returns
$\Omega^{(13)} = 1$.

The Bochner Toeplitz constraint **does** restrict $(a_k, b_k)$ to
genuine Fourier sequences of nonneg measures, but it **does not**
constrain the *higher-order* moments $M[k, k]$, $M[K+k, K+k]$,
$M[k, K+k]$ from above.  These are independent variables in the SDP
lift, with only the PSD lower bounds $M[k,k] \ge a_k^2$,
$M[K+k, K+k] \ge b_k^2$.  The relaxation can still freely inflate
$M[K+k, K+k]$ to absorb the F-side bound on $A_k = a_k^2 - b_k^2$.

The Bochner constraint is **necessary** but not **sufficient** for the
relaxation to give $\Omega > 1$.  We need an additional mechanism
that bounds $M[K+k, K+k]$ from above.

One option (not implemented in v13): use the **dual Bochner
constraint** $[\hat F(j-k)]_{j,k} \succeq 0$ on the lifted F-side
Fourier coefficients $\hat F(k) = M[k,k] - M[K+k,K+k] - 2 i M[k, K+k]$.
This is a separate PSD constraint that mixes the lifted moments.
Whether it helps depends on whether it bounds the diagonal entries
$M[K+k, K+k]$ implicitly.

## Validity (with caveat)

The v13 constraints are **mathematically rigorous**: every admissible
$f$ produces a feasible point with $\Omega^* = \|f * f\|_\infty$.
However, the relaxation is too loose to give $\Omega > 1$.

v14 will add the F-side Bochner constraint on $\hat F$.
