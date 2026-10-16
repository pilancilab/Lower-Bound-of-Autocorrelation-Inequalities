# v7 Proposal: the Dirichlet-kernel Bochner sandwich on $T_f$

## 1. Motivation

v6 imposed
$$0\;\preceq\;T^{f}_{R+1}\;\preceq\;I_{R+1}$$
by applying Bochner to *two* non-negative completions of $f$ to
$[-2,2]$:

* $f$ extended by 0 (yields the lower half $T^f\succeq 0$);
* $\bar f$ defined as $1-f$ on $[-1,1]$ and $1$ outside (yields the
  upper half $I-T^f\succeq 0$).

These are not the only non-negative extensions of $f$-related functions
to $[-2,2]$. Two natural alternatives **also** give non-negative
functions on $[-2,2]$:

* $g_0(x):=(1-f(x))\,\mathbf 1_{[-1,1]}(x)$  —  the *complement* on
  $[-1,1]$ extended by zero outside;
* $f_{\text{ext}}(x):=f(x)\,\mathbf 1_{[-1,1]}(x)+\mathbf 1_{[-2,-1)\cup(1,2]}(x)$
  —  $f$ on $[-1,1]$ extended by $1$ outside.

Both are non-negative on $[-2,2]$, and Bochner-Herglotz gives two new
Toeplitz PSD constraints. Crucially, the Fourier coefficients of
$g_0$ and $f_{\text{ext}}$ at non-zero frequencies are *not* simply
$\pm\hat f$: they pick up an additional **Dirichlet-kernel** term

$$
\mathcal D(j)\;:=\;\frac{\sin(j\pi/2)}{j\pi}\quad (j\ne 0),\qquad \mathcal D(0):=0,
$$
that comes from the $\int_{-1}^{1}e^{-i\pi jx/2}\,dx=4\mathcal D(j)$
contribution of the indicator $\mathbf 1_{[-1,1]}$. The resulting LMIs

$$
\bigl(\tfrac12\bigr)I_{R+1}\;+\;S\;-\;T^{f}_{R+1}\;\succeq\;0,
\qquad
\bigl(\tfrac12\bigr)I_{R+1}\;-\;S\;+\;T^{f}_{R+1}\;\succeq\;0,
$$
where $S_{kl}:=\mathcal D(k-l)$ (with $S_{kk}=0$), are *strictly
sharper* than v6's bounds whenever $\lVert S\rVert_{2}\le \tfrac12$ —
i.e., always for the moderate $R$ we use. These two LMIs together
*sandwich* $T^{f}$ between $S-\tfrac12 I$ and $S+\tfrac12 I$, a window
of width exactly $I$, but **shifted by the Dirichlet kernel**.

## 2. The new constraints

For $R\ge 1$ let the real symmetric $(R+1)\times(R+1)$ Dirichlet
matrix $S$ be
$$
S_{kl}\;=\;\begin{cases}0 & k=l,\\ \dfrac{\sin\bigl((k-l)\pi/2\bigr)}{(k-l)\pi} & k\ne l.\end{cases}
$$

Re-use the v6 expression $T^{f}_{R+1}$ (Hermitian Toeplitz with
diagonal $\tfrac14$ and off-diagonals $(a_{|k-l|}-i\sigma(k-l)b_{|k-l|})/2$).

**Impose**
$$
\boxed{\;\;\tfrac12 I_{R+1}\;+\;S\;-\;T^{f}_{R+1}\;\succeq\;0\quad\text{and}\quad \tfrac12 I_{R+1}\;-\;S\;+\;T^{f}_{R+1}\;\succeq\;0.\;\;}\tag{$\dagger\dagger$}
$$

## 3. Proof

### Step 1. Non-negativity of $g_0$ and $f_{\text{ext}}$ on $[-2,2]$.

On $[-1,1]$: $g_0=1-f\ge 0$ since $f\le 1$; $f_{\text{ext}}=f\ge 0$.
On $[-2,-1)\cup(1,2]$: $g_0\equiv 0\ge 0$; $f_{\text{ext}}\equiv 1\ge 0$.
So both functions are non-negative on $[-2,2]$.

### Step 2. Compute $\hat{g_0}(k)$ and $\hat{f_{\text{ext}}}(k)$.

Both functions agree with $f$ on $[-1,1]$ (modulo replacing $f$ by
$1-f$ for $g_0$) and are constant outside, so the constant-1
contribution gives an indicator-Fourier-coefficient
$\hat{\mathbf 1}_{[-1,1]}(k)$ at the appropriate sign.

A direct computation: $\int_{-1}^{1}e^{-i\pi kx/2}\,dx$
$=\bigl[-\tfrac{2}{i\pi k}e^{-i\pi kx/2}\bigr]_{-1}^{1}
=\tfrac{4\sin(k\pi/2)}{k\pi}=4\,\mathcal D(k)$ for $k\ne 0$,
and $=2$ for $k=0$.

**$\hat{g_0}(0)$**: $\hat{g_0}(0)=\tfrac14\int_{-1}^{1}(1-f)\,dx=\tfrac14\cdot 1=\tfrac14$.

**$\hat{g_0}(k)$ for $k\ne 0$**:
$$
\hat{g_0}(k)=\tfrac14\!\int_{-1}^{1}\!(1-f)e^{-i\pi kx/2}dx
=\tfrac14\bigl[4\mathcal D(k)-4\hat f(k)\bigr]
=\mathcal D(k)-\hat f(k).
$$

**$\hat{f_{\text{ext}}}(0)$**: $\hat{f_{\text{ext}}}(0)=\tfrac14\bigl[\int_{-1}^{1}\!f+\int_{[-2,-1)\cup(1,2]}\!1\bigr]=\tfrac14\bigl[1+2\bigr]=\tfrac34$.

**$\hat{f_{\text{ext}}}(k)$ for $k\ne 0$**: writing $f_{\text{ext}}=f\cdot\mathbf 1_{[-1,1]}+\mathbf 1_{[-2,2]}-\mathbf 1_{[-1,1]}$ and noting
$\hat{\mathbf 1_{[-2,2]}}(k)=0$ for $k\ne 0$ integer,
$$
\hat{f_{\text{ext}}}(k)=\hat f(k)+0-\mathcal D(k)=\hat f(k)-\mathcal D(k).
$$

### Step 3. Identify the new Toeplitz matrices.

For $T^{g_0}_{R+1}=[\hat{g_0}(k-l)]_{0\le k,l\le R}$:

* diagonal $\tfrac14$;
* off-diagonal $\mathcal D(k-l)-\hat f(k-l)$.

Compared to $T^{f}_{R+1}$, this is

$$
T^{g_0}_{R+1}\;=\;\tfrac12 I_{R+1}\;+\;S\;-\;T^{f}_{R+1},
$$

since the diagonal becomes $\tfrac14=\tfrac12+0-\tfrac14$ (using
$S_{kk}=0$), and the off-diagonals become $0+\mathcal D(k-l)-\hat f(k-l)$.

For $T^{f_{\text{ext}}}_{R+1}=[\hat{f_{\text{ext}}}(k-l)]$:

* diagonal $\tfrac34$;
* off-diagonal $\hat f(k-l)-\mathcal D(k-l)$.

A parallel calculation gives

$$
T^{f_{\text{ext}}}_{R+1}\;=\;\tfrac12 I_{R+1}\;-\;S\;+\;T^{f}_{R+1}.
$$

### Step 4. Bochner.

By the same trigonometric polynomial argument as v4/v6, since
$g_0\ge0$ and $f_{\text{ext}}\ge 0$ on $[-2,2]$,
$$
T^{g_0}_{R+1}\;\succeq\;0\quad\text{and}\quad T^{f_{\text{ext}}}_{R+1}\;\succeq\;0,
$$
which by Step 3 are exactly the two LMIs in $(\dagger\dagger)$.

### Step 5. Admissible $\Rightarrow$ feasible.

At any admissible $f$, the program-point assignment gives
$a_k^{\text{prog}}=a_k^{\text{true}}$ and $b_k^{\text{prog}}=b_k^{\text{true}}$
(as in v6), so $T^{f}_{R+1}$ equals the true Toeplitz $[\hat f(k-l)]$.
At that point $g_0$ and $f_{\text{ext}}$ have the form described above
and are nonneg on $[-2,2]$, so by Step 4 $(\dagger\dagger)$ holds. All
v0–v6 constraints continue to hold as before. ✓

### Step 6. Convexity.

$T^{f}_{R+1}$ is affine in the program variables $(a_k,b_k)$ (v6
proof). $S$ is a *constant* matrix. So both
$\tfrac12 I+S-T^{f}$ and $\tfrac12 I-S+T^{f}$ are affine Hermitian
matrix-valued functions, and the constraints $\succeq 0$ are LMIs —
convex. Implementation in the real form (the same
$(2R+2)\times(2R+2)$ block trick as v6) is standard. ✓

## 4. Why this is genuinely new

* No constraint in v0–v6 involves the Dirichlet kernel
  $\mathcal D(j)=\sin(j\pi/2)/(j\pi)$. The matrix $S$ is the
  truncated Toeplitz of the indicator $\mathbf 1_{[-1,1]}$ on
  $[-2,2]$ — a *non-trivial constant* that encodes how the support
  of $f$ sits inside the period-4 interval $[-2,2]$.
* The new sandwich $S-\tfrac12 I\preceq T^{f}\preceq S+\tfrac12 I$
  shifts the v6 sandwich $0\preceq T^{f}\preceq I$ *off-center by $S$*.
  In PSD-cone terms this is strictly different from v6 — neither one
  implies the other in general.
* In particular, at the test regime where $\hat f(1)\approx p/2\approx 0.19$
  (real, since $d_1\approx 0$), the (0,1) entry of the new "upper"
  matrix is $\mathcal D(1)+\tfrac12-\hat f(1)=\tfrac{1}{\pi}-\hat f(1)\approx 0.318-0.19=0.128$
  vs v6's $\tfrac12-\hat f(1)\approx 0.31$.  Strictly tighter.

## 5. Strictness and combination with v6

v6 + v7 together impose
$$
\max(0,\,S-\tfrac12 I)\;\preceq\;T^{f}_{R+1}\;\preceq\;\min(I,\,S+\tfrac12 I)
$$
(in the PSD sense, with max/min interpreted as the natural
intersection of the relevant cones). For moderate $R$ where
$\lVert S\rVert_{2}\le\tfrac12$, $S-\tfrac12 I\preceq 0$ and
$S+\tfrac12 I\preceq I$, so both new LMIs are strictly tighter than
the respective v6 LMIs whenever $S\ne 0$ — which is always (the
Dirichlet kernel has non-zero entries at odd offsets).

## 6. Implementation note

The constant matrix $S$ is computed once at program-construction time
as a numpy array of shape $(R+1,R+1)$ with $S_{kl}=\sin((k-l)\pi/2)/((k-l)\pi)$
for $k\ne l$ and $S_{kk}=0$. The real-form embedding of $S$ in the
$(2R+2)\times(2R+2)$ block reformulation is the block-diagonal
$\bigl(\begin{smallmatrix}S & 0\\ 0 & S\end{smallmatrix}\bigr)$
(since $S$ is real and symmetric, its real form is itself,
duplicated). The two new LMIs are then
$$
\tfrac12 I_{2R+2}+\mathrm{blkdiag}(S,S)-S^{f}_{\text{real}}\succeq 0,
\qquad
\tfrac12 I_{2R+2}-\mathrm{blkdiag}(S,S)+S^{f}_{\text{real}}\succeq 0,
$$
each $(2R+2)\times(2R+2)$ real symmetric PSD constraints — same size
as v6's LMIs. No new optimisation variables are needed.
