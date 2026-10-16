# v6 Proposal: Bochner sandwich on $f$ itself

## 1. Motivation

v4 and v5 applied Bochner–Herglotz to $M$ and to $\Omega-M$. Both
operators live on the *output* side of the problem. The **input**
side — the function $f$ itself — has so far been touched only through
its $L^2$ energy ($\|f\|_2^2\le 1$, i.e. (5.11) of White) and through
the convolution identity that links $\hat f$ to $\hat M$.

But the original problem specification says
$$f:[-1,1]\to[0,1],\qquad g=1-f\quad\text{on }[-1,1],$$
which means *both*
$$f(x)\ge 0\quad\text{and}\quad 1-f(x)\ge 0\quad\text{on }[-1,1].$$
With the extension-by-zero convention used throughout White's paper,
$f$ and $g$ are both non-negative on $[-2,2]$ as well: outside $[-1,1]$
each is identically $0$. Therefore *both* $f\ge 0$ and the function
$$\bar f(x):=
\begin{cases}1-f(x) & x\in[-1,1],\\ 1 & x\in[-2,-1)\cup(1,2]\end{cases}$$
are non-negative on $[-2,2]$. Bochner on each yields a Toeplitz PSD
constraint. Together they sandwich the Toeplitz of $\hat f$ between
$0$ and $I$ — analogously to v5's sandwich on $T_M$, but now on the
$f$-side of the problem. No version up to v5 imposes these.

## 2. The new constraints

For each $k\in\{0,1,\dots,R\}$ let $\hat f(k)$ denote the Fourier
coefficient of $f$ on $[-2,2]$. By White's (3.2)
$\hat f(k)=(a_k-i b_k)/2$ for $k\ge 1$, and
$\hat f(0)=\tfrac14\!\int_{-2}^{2}\!f\,dx=\tfrac14$ (using
$\int_{-1}^{1}\!f=1$ and $f\equiv 0$ off $[-1,1]$). Define the
Hermitian Toeplitz expression
$$
T_{R+1}^{f}\;:=\;\bigl[\hat f(k-l)\bigr]_{0\le k,l\le R}\quad\text{with}\quad T_{kk}^{f}=\tfrac14,\quad T_{kl}^{f}=\tfrac{a_{|k-l|}-i\sigma(k-l)b_{|k-l|}}{2}.
$$

Impose
$$
\boxed{\;T_{R+1}^{f}\;\succeq\;0\quad\text{and}\quad I_{R+1}\;-\;T_{R+1}^{f}\;\succeq\;0.\;} \tag{$\dagger$}
$$

## 3. Proof

### Step 1. $f\ge 0$ and $\bar f\ge 0$ on $[-2,2]$.

Direct from the problem statement plus the extension-by-zero
convention used by White. $f$ is non-negative on $[-1,1]$ and $0$
outside; $\bar f$ equals $g=1-f\ge 0$ on $[-1,1]$ and $1\ge 0$ outside.

### Step 2. Compute $\hat f(0)$ and $\hat{\bar f}(0)$.

$\hat f(0)=\tfrac14\!\int_{-2}^{2}\!f\,dx=\tfrac14\cdot 1=\tfrac14$.

$\hat{\bar f}(0)=\tfrac14\!\int_{-2}^{2}\!\bar f\,dx
=\tfrac14\!\left[\int_{-1}^{1}\!(1-f)\,dx + \int_{[-2,-1)\cup(1,2]}\!1\,dx\right]
=\tfrac14\!\left[(2-1) + 2\right]=\tfrac{3}{4}$.

### Step 3. Compute $\hat{\bar f}(k)$ for $k\ne 0$.

For $k\ne 0$ integer, $\int_{-2}^{2}\!e^{-i\pi kx/2}\,dx=0$. So
$$
\hat{\bar f}(k)
=\tfrac14\!\int_{-2}^{2}\!\bar f(x)e^{-i\pi kx/2}\,dx
=\tfrac14\!\int_{-2}^{2}\!e^{-i\pi kx/2}\,dx
\;-\;\tfrac14\!\int_{-2}^{2}\!f(x)e^{-i\pi kx/2}\,dx
=0-\hat f(k)
=-\hat f(k).
$$
Therefore the Toeplitz matrix $T^{\bar f}$ has entries
$$
T^{\bar f}_{kl}=\hat{\bar f}(k-l)=\begin{cases}3/4 & k=l,\\ -\hat f(k-l) & k\ne l,\end{cases}
$$
which is exactly $I_{R+1}-T^{f}_{R+1}$ (the off-diagonals contribute
$-\hat f$ and the diagonal goes from $1/4$ to $3/4=1-1/4$).

### Step 4. Bochner.

For any $u_0,\dots,u_R\in\mathbb C$, the trigonometric polynomial
$p(x)=|\sum u_k e^{i\pi kx/2}|^{2}$ is non-negative on $\mathbb R$.
Then
$$0\le\int_{-2}^{2}p(x)f(x)\,dx
=4\sum_{k,l}u_k\overline{u_l}\,\overline{\hat f(k-l)},$$
giving $T^{f}_{R+1}\succeq 0$ (as in v4's Step 1, with $f$ in place of
$M$). By exactly the same argument with $\bar f$ in place of $f$,
$T^{\bar f}_{R+1}=I-T^{f}_{R+1}\succeq 0$. So both halves of $(\dagger)$
hold at every admissible $f$.

### Step 5. Convexity.

The entries of $T^{f}_{R+1}$ are *affine* in the program variables
$a_k,b_k$ — the diagonal is the constant $1/4$, and the off-diagonals
are $(a_{|k-l|}\pm i b_{|k-l|})/2$, both linear in $(a_{|k-l|},b_{|k-l|})$
which are themselves program variables (already constrained by Lemma 3
to be the truncated/exact Fourier coefficients of $f$ on $[-2,2]$).
So $T^{f}_{R+1}$ is an affine Hermitian matrix-valued function, and
both $T^{f}\succeq 0$ and $I-T^{f}\succeq 0$ are linear matrix
inequalities — convex. The standard real reformulation
$\bigl(\begin{smallmatrix}\Re T & -\Im T\\ \Im T & \Re T\end{smallmatrix}\bigr)\succeq 0$
is used for the implementation, giving a $(2R+2)\times(2R+2)$ real
symmetric PSD constraint of the same size as v4's. ✓

### Step 6. Admissible $\Rightarrow$ feasible.

At any admissible $f$, with the program point assigned via
$c_k,d_k$ and the resulting $a_k,b_k$ from Lemma 3 (3.3)/(3.4)
(including the Lemma 4 tail $\varepsilon_k,\delta_k$ choice that makes
$a_k^{\text{prog}}=a_k^{\text{true}}$ and similarly for $b_k$),
the Toeplitz expression $T^{f}_{R+1}$ in the program becomes the *true*
Toeplitz of $\hat f$. By Step 4 it is PSD and $I-T^{f}\succeq 0$.
All v0–v5 constraints continue to hold as before. So every admissible
$f$ extends to a v6-feasible point and the program optimum continues
to lower-bound $\inf_f\lVert M\rVert_\infty$.

## 4. Why this is genuinely new

* All previous v0–v5 constraints on the *input* function $f$ are
  *energy* constraints in disguise: (5.10) is an $L^\infty$-on-Fourier-
  coefficients bound, (5.11) is the Parseval $L^2$ bound, and Lemmas 2
  and 3 are *exact identities* (no inequality). No previous constraint
  imposes a positivity / Bochner-type *structural* property on
  $\hat f$. v6 fills that gap.
* The Toeplitz PSD on $\hat f$ encodes *correlation inequalities* among
  the Fourier coefficients $\hat f(k)$ at different $k$ — strictly more
  than the per-coefficient or per-pair-energy bounds of v0–v5.
* The dual half $I-T^{f}\succeq 0$ is the Bochner-side translation of
  the upper bound $f\le 1$ — completely missing from v0–v5 except in
  the very weak forms $|c_k|\le 2/\pi$ and the implicit
  $\sum(c_k^{2}+d_k^{2})\le 1/2$ (Parseval applied to $f^2\le f$).

## 5. Strictness

v6 ⊃ v5: every v5-feasible point assigns $a_k,b_k$ consistent with
some $\hat f$, and adding $(\dagger)$ removes those points where the
resulting Toeplitz fails PSD or exceeds $I$. So
$\Omega_{\text{v6}}^{\star}\ge\Omega_{\text{v5}}^{\star}$, with strict
gain unless v5 happens to already produce $a,b$ corresponding to a
genuinely positive density bounded by $1$ — which is the *defining
condition* on $f$ that v0–v5 only enforced through energy shadows.

## 6. Implementation note

The Toeplitz on $\hat f$ uses the program variables $a_k,b_k$
**directly** — no further lifting needed. The two new LMIs are the
real-form PSD constraints
$$
\begin{pmatrix}\Re T^{f} & -\Im T^{f}\\ \Im T^{f} & \Re T^{f}\end{pmatrix}\succeq 0
\qquad\text{and}\qquad
\begin{pmatrix}\Re(I-T^{f}) & -\Im(I-T^{f})\\ \Im(I-T^{f}) & \Re(I-T^{f})\end{pmatrix}\succeq 0,
$$
both of size $(2R+2)\times(2R+2)$. With $R=10$ each is $22\times 22$
— same scale as the v4/v5 SDPs and well within MOSEK/SCS capabilities.
