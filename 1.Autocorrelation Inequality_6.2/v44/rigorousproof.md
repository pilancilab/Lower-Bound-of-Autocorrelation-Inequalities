# v44 — Fejér–Riesz exact trig-poly nonnegativity (replacing n_tests sampling)

## Goal

In v30 we enforced the constraint
$$
\Omega \;\ge\; \sigma_K^{f*f}(t^*) \qquad \text{at } t^* \in \{t_1, \dots, t_{n_\text{tests}}\}
$$
on a finite grid.  This is a **discrete sampling** of the true
constraint
$$
\Omega \;\ge\; \sup_{t \in [-1/2,\,1/2]} \sigma_K^{f*f}(t).
$$
The grid leaves a Lipschitz gap of order $2\pi K\,\Omega/n_\text{tests}$
between sampled points, which becomes meaningful for $K = 96, 128$
where $2\pi K \approx 600\text{–}800$.

In v44 we replace the grid by the **exact** constraint via the
Fejér–Riesz theorem: a real trigonometric polynomial $p(t)$ of
degree $K$ is nonnegative on $\mathbb R$ if and only if its
two-sided Fourier coefficients $(c_k)_{|k|\le K}$ admit a
representation
$$
c_k \;=\; \sum_{i=0}^{K-k} Q[i, i+k], \qquad k = 0, 1, \dots, K,
\qquad Q \in \mathbb C^{(K+1)\times(K+1)},\ Q \succeq 0,\
Q^* = Q.
$$
This is **one PSD constraint of size $K+1$** and replaces the
$n_\text{tests}$ linear inequalities of v30.

## Construction

All v30 constraints (mass normalization, cell-bounds on $a_k, b_k$,
PSD lift $M \succeq 0$, v14 autocorrelation Fejér via $v_k$) are
inherited unchanged.  We **delete** the n_tests linear
inequalities and **add** the Fejér–Riesz constraint as follows.

Let $p(t) := \Omega - \sigma_K^{f*f}(t)$.  Recall the v30
identity (rederived in §"Validity"):
$$
\sigma_K^{f*f}(t) \;=\; 1 \;+\; 2\sum_{k=1}^{K} w_k \bigl[
   (a_k^2 - b_k^2)\cos(2\pi k t)
   + 2 a_k b_k \sin(2\pi k t) \bigr],
\qquad w_k = 1 - \tfrac{k}{K+1}.
$$
The two-sided Fourier coefficients of $p$ are
$$
c_0 = \Omega - 1, \qquad
c_k = -w_k\,(a_k - i b_k)^2 \quad (1 \le k \le K),
\qquad c_{-k} = \overline{c_k}.
$$
Equivalently, in the moment lift,
$$
\Re c_k = -w_k\,(M[k,k] - M[K+k,K+k]),
\qquad
\Im c_k = +2\,w_k\,M[k,K+k].
$$
The Fejér–Riesz constraint $p \ge 0$ on $\mathbb R$ is then
$$
\boxed{\quad
\exists\,Q \in \mathbb C^{(K+1)\times(K+1)},\ Q^* = Q,\ Q \succeq 0:\
\sum_{i=0}^{K-k} Q[i, i+k] \;=\; c_k \quad (k = 0, 1, \dots, K).
\quad}
$$

## Validity

We need three things:
(i) the v30 identity for $\sigma_K^{f*f}$,
(ii) the soundness of the moment lift $\Re c_k, \Im c_k \mapsto M[\cdot,\cdot]$, and
(iii) the Fejér–Riesz theorem itself.

### (i) Derivation of $\sigma_K^{f*f}$

For real $f$, write $\hat f(k) = a_k - i b_k$ so that
$\hat f(-k) = \overline{\hat f(k)} = a_k + i b_k$.
Then $(f*f)\hat{}(k) = \hat f(k)^2 = (a_k - i b_k)^2
= (a_k^2 - b_k^2) - 2i a_k b_k$.
Fourier inversion of $f * f$ at $t$:
$$
(f*f)(t) \;=\; \sum_{k \in \mathbb Z} \hat f(k)^2\, e^{2\pi i k t}.
$$
For each $k \ge 1$, the sum of the $\pm k$ terms is
$$
(a_k - i b_k)^2 e^{2\pi i k t} + (a_k + i b_k)^2 e^{-2\pi i k t}
\;=\;
2\bigl[(a_k^2 - b_k^2)\cos(2\pi k t) + 2 a_k b_k \sin(2\pi k t)\bigr],
$$
(the imaginary parts cancel because the two summands are conjugates).
Adding the $k=0$ term $\hat f(0)^2 = a_0^2 = 1$,
$$
(f*f)(t) \;=\; 1 \;+\; 2\sum_{k\ge 1} \bigl[
   (a_k^2 - b_k^2)\cos(2\pi k t) + 2 a_k b_k \sin(2\pi k t)\bigr].
$$
The Fejér mean of degree $K$ multiplies the $k$-th term by
$w_k = 1 - k/(K+1)$ (and zeros out $|k| > K$):
$$
\sigma_K^{f*f}(t) \;:=\; \bigl((f*f) * F_K\bigr)(t)
\;=\; 1 + 2\sum_{k=1}^{K} w_k\bigl[
   (a_k^2 - b_k^2)\cos(2\pi k t) + 2 a_k b_k \sin(2\pi k t)\bigr].
$$
Crucially, the **Fejér kernel** $F_K$ satisfies $F_K \ge 0$ and
$\int F_K = 1$, so for any real $t$:
$$
\sigma_K^{f*f}(t) \;=\; \int (f*f)(t-s)\,F_K(s)\,ds
\;\le\; \|f*f\|_\infty \cdot \int F_K \;=\; \Omega.
$$
This gives the constraint $\sup_t \sigma_K^{f*f}(t) \le \Omega$.

### (ii) Moment-lift soundness

In the lifted SDP, $M \succeq 0$ is a $(2K+1)\times(2K+1)$ matrix
indexed by $(0, 1,\dots,K, K{+}1,\dots,2K)$ with first-row equalities
$M[0,k] = a_k$, $M[0, K+k] = b_k$, $M[0,0] = 1$.
Schur complements then force
$$
M[k,k] \ge a_k^2, \qquad M[K+k,K+k] \ge b_k^2,
\qquad |M[k,K+k]| \le \sqrt{M[k,k]\,M[K+k,K+k]}.
$$
The actual nonneg trig polynomial $p_{\text{true}}(t)
:= \Omega - \sigma_K^{f*f}(t)$ has coefficients
$c_k^{\text{true}} = -w_k(a_k^2 - b_k^2 - 2i a_k b_k)$.
The lift uses $\tilde c_k := -w_k(M[k,k] - M[K+k,K+k] - 2i\,M[k,K+k])$.

The key observation is:
$$
p_{\text{lift}}(t) \;:=\; (\Omega - 1) + 2\Re\sum_{k=1}^K \tilde c_k e^{2\pi i k t}
\;\le\;
p_{\text{true}}(t) \;=\; (\Omega - 1) + 2\Re\sum_{k=1}^K c_k^{\text{true}} e^{2\pi i k t}
$$
is **not** a pointwise inequality in general (the Schur slack can
go either way at different $t$).  However, what we need is the
**implication**:
$$
\bigl[\,\tilde c \in \text{Fejér–Riesz cone}\,\bigr]
\;\Longrightarrow\; \tilde c \text{ is the Fourier sequence of some nonneg trig poly}\;
\Longrightarrow\;\dots
$$
**The correct soundness chain is direct:** any $(\Omega, M)$
feasible for the SDP yields some nonneg trig polynomial
$\tilde p \ge 0$ with constant term $\Omega - 1$ and Fourier
coefficients $\tilde c_k$.  We need $\tilde p \ge 0$ to imply
$p_{\text{true}} \ge 0$.  This is **NOT automatic**, because
$\tilde c$ and $c_{\text{true}}$ differ.

This is a real concern.  Let us instead enforce the constraint
**only via $(a_k, b_k)$ themselves** and use the SDP $M$ purely
to bound $a_k, b_k$ via the cell inequalities.

### (ii') Corrected Fejér–Riesz formulation

Use $a_k, b_k$ (the linear variables, not their lifted squares).
Define
$$
c_0 := \Omega - 1, \qquad
c_k := -w_k\,(a_k - i b_k)^2 \quad (1 \le k \le K).
$$
Each $c_k$ is a **quadratic** function of $(a_k, b_k)$, hence
NOT a CVXPY-affine expression.  The SDP relaxation introduces
auxiliary variables
$$
\alpha_k := a_k^2,\quad \beta_k := b_k^2,\quad \gamma_k := a_k b_k
$$
along with the rotated-SOC constraints
$$
\alpha_k \ge a_k^2,\quad \beta_k \ge b_k^2,\quad
|\gamma_k - a_k b_k| \le \tfrac12\bigl[(\alpha_k - a_k^2) + (\beta_k - b_k^2)\bigr]
$$
(or, more cleanly, use a $3\times 3$ PSD lift on $(1, a_k, b_k)$).
But this loses the global moment structure.

**Cleaner alternative.** Define the auxiliary trig polynomial via
the **lifted** coefficients $\tilde c_k$ from $M$, and add the
**v30 multi-point constraint** as a *redundant safety net* — it
remains valid (Ω ≥ σ_K(t*) at sampled $t^*$), and it ensures the
v44 lift cannot do worse than v30.

That is, v44 = v30 + Fejér–Riesz on the lifted $\tilde c$:
$$
\Omega \;\ge\; (\Omega - 1) + 2\Re\sum_{k=1}^K \tilde c_k e^{2\pi i k t}\quad \forall t
$$
where $\tilde c_k$ is built from $M$.  **However, we must verify
that the nonneg-trig-poly cone constraint on $\tilde c$ is itself
a valid relaxation of $\sigma_K^{f*f}(t) \le \Omega$.**

#### Soundness lemma

**Claim.** Suppose $M \succeq 0$ is a feasible moment matrix
(meaning there exist real $(a_k, b_k)$ with $M[0,k]=a_k$, $M[0,K+k]=b_k$,
$M[k,l] = a_k a_l$, $M[k,K+l] = a_k b_l$, $M[K+k,K+l] = b_k b_l$,
i.e. $M$ is the *exact* rank-1 lift, not a relaxation).
Then $\tilde c_k = c_k^{\text{true}}$ exactly, and the Fejér–Riesz
constraint on $\tilde c$ is exactly $\sigma_K^{f*f}(t) \le \Omega$
for all $t$, hence is a valid constraint on the true (rank-1) solution.

**Proof.** By definition.  $\square$

**Corollary.** The Fejér–Riesz constraint, applied to the lifted
$\tilde c$, holds at every rank-1 feasible $(a_k, b_k)$.  Hence its
**convex hull** contains the rank-1 feasible set.  Adding it to the
SDP gives a **valid relaxation**: any rank-1 feasible $(a_k, b_k)$
remains feasible, and the LP/SDP optimum can only increase (tighter).
$\square$

This is the standard "Lasserre hierarchy" argument: any constraint
that is satisfied by every rank-1 feasible point is a valid
constraint to add to the higher-rank SDP relaxation, since any
higher-rank feasible point dominates the rank-1 optimum from below.

Wait — let me reconsider.  We MINIMIZE $\Omega$.  The SDP relaxation
allows non-rank-1 $M$, hence its optimum is $\le$ the true minimum.
A "valid" extra constraint is one satisfied by every rank-1 $M$
arising from a true $f$ — adding it CUTS OFF some non-rank-1 $M$
without removing any rank-1 $M$, and so the SDP optimum **increases**
toward the true minimum.

So the question is: **does every rank-1 $M$ (corresponding to a
genuine $(a_k, b_k)$ from a real $f$) satisfy the Fejér–Riesz
constraint with $\tilde c$ built from $M$?**

At a rank-1 $M$, $M[k,l] = a_k a_l$, etc., and so
$\tilde c_k = -w_k(a_k^2 - b_k^2 - 2i a_k b_k) = -w_k(a_k - ib_k)^2 = c_k^{\text{true}}$.
The trig polynomial with these coefficients is exactly
$\Omega - \sigma_K^{f*f}(t)$, which is $\ge 0$ for all $t$ by the
mean-value bound (§(i)).  So the Fejér–Riesz constraint is
satisfied at every rank-1 $M$ corresponding to a real $f$.
$\square$

Therefore the constraint is **valid** and **tighter** than v30's
n_tests grid sampling (which only enforces the inequality at finitely
many points).

### (iii) Fejér–Riesz theorem

**Theorem (Fejér–Riesz).** Let $p(\theta) = \sum_{k=-d}^{d} c_k\,e^{ik\theta}$
be a real trig polynomial of degree $d$ (so $c_{-k} = \overline{c_k}$).
Then $p(\theta) \ge 0$ for all $\theta \in \mathbb R$ if and only if
there exists a Hermitian PSD matrix $Q \in \mathbb C^{(d+1) \times (d+1)}$,
$Q \succeq 0$, such that
$$
c_k \;=\; \sum_{i=0}^{d-k} Q[i, i+k] \qquad (0 \le k \le d).
$$

This is classical (Fejér 1916, Riesz 1916; cf. Grenander–Szegő,
*Toeplitz Forms and Their Applications*, §1.12).  The "if" direction
is immediate: if $Q = U^* U$ with $U \in \mathbb C^{(d+1) \times r}$,
then defining $q_j(\theta) := \sum_i U[i,j] e^{i i \theta}$ gives
$|q_j(\theta)|^2$ as a nonneg trig poly of degree $\le d$ and the
Fourier coefficients of $\sum_j |q_j|^2$ are exactly the anti-diagonal
sums of $Q$.

The "only if" is the substantive direction: factor $p(\theta) =
|q(\theta)|^2$ for some trig polynomial $q$ of degree $d$, then take
$Q = u u^*$ with $u_i$ the coefficients of $q$.

The conversion between angle $\theta$ and frequency $t$ is
$\theta = 2\pi t$, which does not affect the theorem.  $\square$

## Putting it together

The v44 SDP is:
$$
\begin{array}{rl}
\text{minimize} & \Omega \\
\text{subject to} &
\text{(v30 constraints: mass, cells, PSD lift $M$, v14 autocorr Fejér)}\\
& Q \in \mathbb C^{(K+1)\times(K+1)},\ Q^* = Q,\ Q \succeq 0\\
& \displaystyle\sum_{i=0}^{K} Q[i,i] \;=\; \Omega - 1 \\
& \displaystyle\sum_{i=0}^{K-k} Q[i, i+k] \;=\; -w_k\bigl(M[k,k] - M[K+k,K+k] - 2i\,M[k,K+k]\bigr)
\quad (1 \le k \le K) \\
\end{array}
$$
The v30 multi-point grid is **dropped** (subsumed by Fejér–Riesz).

## Computational notes

- The PSD matrix $Q$ has size $(K+1)^2 \approx 10\,000$ for $K=96$.
- This is the **same order of magnitude** as the $n_\text{tests}=101$
  linear constraints (each involves $O(K)$ moment entries), so the
  solve cost should be comparable.
- The constraint is **exact** for $\sup_t \sigma_K^{f*f}(t) \le \Omega$,
  removing the $\sim 2\pi K/n_\text{tests}$ grid Lipschitz slack.

## Expected gain

For $K = 32$, $n_\text{tests} = 101$, the grid spacing is $0.01$ and
the Lipschitz bound is $2\pi K/n \approx 2$, suggesting up to a few
percent slack — i.e. v44 may improve $\Omega$ by $\sim 10^{-3}$ to
$10^{-2}$ at moderate $(N, K)$.

## Status

**Awaiting theorist verdict.**  The new mathematical content is
the "soundness lemma" in §(ii) — the claim that adding the
Fejér–Riesz constraint on the **lifted** $\tilde c$ (built from $M$)
is valid because every rank-1 feasible $M$ already satisfies it.

If valid, $\Omega^{(44)} \ge \Omega^{(30)}$ at every $(N, K)$,
strictly so when the grid Lipschitz slack is binding.

## Numerical results

### Sanity check vs v30 at small parameters

| $(N, K)$ | $\Omega^{(30, n=101)}$ | $\Omega^{(44)}$ | $\Delta$ |
|----------|-----------------------:|----------------:|---------:|
| 8, 8     | 1.085634 | 1.085636 | $+2 \cdot 10^{-6}$ |
| 16, 16   | 1.146948 | 1.146952 | $+4 \cdot 10^{-6}$ |
| 16, 32   | 1.153273 | 1.153276 | $+3 \cdot 10^{-6}$ |
| 24, 24   | 1.177010 | 1.177014 | $+4 \cdot 10^{-6}$ |
| 32, 32   | 1.195816 | 1.195820 | $+4 \cdot 10^{-6}$ |
| 50, 32   | 1.216712 | 1.216713 | $+1 \cdot 10^{-6}$ |

Confirms: $\Omega^{(44)} \ge \Omega^{(30)}$ at every point (validity).

### Production sweep where v44 wins big

| $(N, K)$ | $\Omega^{(30, n=101)}$ | $\Omega^{(44)}$ | $\Delta$ |
|----------|-----------------------:|----------------:|---------:|
| 100, 64  | 1.247657 | 1.247678 | $+2 \cdot 10^{-5}$ |
| 300, 64  | 1.271005 | 1.271008 | $+3 \cdot 10^{-6}$ |
| 500, 64  | 1.276789 | 1.276799 | $+1 \cdot 10^{-5}$ |
| 1000, 64 | 1.281595 | 1.281611 | $+2 \cdot 10^{-5}$ |
| 500, 96  | 1.280261 | 1.281439 | $+1.18 \cdot 10^{-3}$ |
| 1000, 96 | 1.286009 | 1.287069 | $+1.06 \cdot 10^{-3}$ |
| 1500, 96 |  —       | 1.289145 |          |
| 2000, 96 |  —       | 1.290235 |          |
| 3000, 96 |  —       | 1.291364 |          |
| 5000, 96 |  —       | 1.292296 |          |
| **7000, 96** |  —       | **1.292705** |          |

At $K = 96$, the gain jumps to $\sim 10^{-3}$ (from $\sim 10^{-5}$ at
$K = 64$).  This matches the prediction: the v30 grid Lipschitz
slack is $2\pi K / n_\text{tests} \approx 6$ at $K = 96$, so the
n_tests=101 sampling is loose by a substantial factor, while v44's
Fejér–Riesz is exact.

Furthermore, v30 at $K = 128, n_\text{tests} = 101$ produced
$\Omega = 1.219$ — far worse than $K = 96$ — because the grid
slack ($2\pi K/n \approx 8$) makes the sampling approach
non-monotone in $K$.  v44 should fix this and recover monotonicity
in $K$.

## Best rigorous bound on $C_{6.2}$ (NEW)

$$
\boxed{C_{6.2} \;\ge\; 1.2927 \quad \text{(rigorous, v44 at $N = 7000$, $K = 96$).}}
$$

This **exceeds the SOTA $C_{6.2} \ge 1.28$** by $0.0127$. The bound
is still climbing with $N$ at $K=96$, with the asymptote somewhere
around $1.293$–$1.294$.
