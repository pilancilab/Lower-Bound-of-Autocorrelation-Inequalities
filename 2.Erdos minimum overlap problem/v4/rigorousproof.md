# v4 Proposal: Bochner–Herglotz Toeplitz PSD constraint

## 1. Motivation

The v1–v3 chain extracted everything the *energy* of $M$ has to say:
$\int M^{2}\le \Omega$ plus Lemmas 3, 4, 5 gave us a sum-of-squares
bound on the Fourier coefficients of $M$, and v3 closed it by lifting
the odd-index cosine block into a convex variable.

But energy is only one shadow of the deeper fact

$$M(x)\;\ge\;0\quad\text{everywhere on }[-2,2].$$

A positive measure on a finite interval has Fourier coefficients
$\{\hat M(k)\}_{k\in\mathbb Z}$ forming a **positive-definite sequence**
(Bochner–Herglotz). Equivalently, every finite Hermitian Toeplitz
matrix
$$
T_{R+1}\;:=\;\bigl[\hat M(k-l)\bigr]_{0\le k,l\le R}
$$
is positive semi-definite. Energy gives only $\|\hat M\|_{\ell^{2}}^{2}$
bounds — $T_{R+1}\succeq 0$ is a *strictly stronger* family of
constraints, including all the principal-minor inequalities (which are
multilinear in $\hat M(k)$ and capture *correlations* between Fourier
coefficients of $M$, not just their magnitudes).

## 2. The new constraint

Extend the v3 *lifting* to all $k\in\{1,2,\dots,R\}$:

For every $k\in\{1,\dots,R\}$, introduce a free variable
$\tilde A_{k}\in\mathbb R$ with

$$
F_{\cos}^{\text{LB}}[k]\;\le\;\tilde A_{k}\;\le\;F_{\cos}^{\text{UB}}[k], \tag{A1}
$$
$$
\tilde A_{k}\;\le\;\frac{4\sin(k\pi/2)}{k\pi}\,a_{k}\;-\;2\bigl(a_{k}^{2}+b_{k}^{2}\bigr). \tag{A2}
$$
For odd $k=2r-1$ this is exactly the v3 lifting; for even $k=2m$, (A2)
becomes $\tilde A_{2m}\le -2(a_{2m}^{2}+b_{2m}^{2})=-\tfrac12(c_{m}^{2}+d_{m}^{2})$
(since $\sin(m\pi)=0$ and $a_{2m}=c_{m}/2,\ b_{2m}=d_{m}/2$).

Now define the $(R+1)\times(R+1)$ **Hermitian Toeplitz expression**
$$
T_{k,l}\;:=\;
\begin{cases}
\dfrac14 & k=l,\\[1mm]
\dfrac{\tilde A_{|k-l|}-i\,\sigma(k-l)\,B_{|k-l|}}{2} & k\neq l,
\end{cases}
$$
where $\sigma(\cdot)\in\{\pm 1\}$ is the sign of its argument, and
$B_{j}$ is the program expression for the sine Fourier coefficient of
$M$:
$$
B_{j}\;:=\;\begin{cases}
-\dfrac{4(-1)^{r+1}}{j\pi}\,b_{j}^{\text{prog}} & j=2r-1\text{ odd},\\
0 & j\text{ even}
\end{cases}
$$
(matching Lemma 3 (3.6); the program variable $b_{j}^{\text{prog}}$ for
odd $j$ already encodes the true sine Fourier coefficient of $f$ on
$[-2,2]$, as established in the v2 verdict).

**Impose**
$$
\boxed{\;T_{R+1}\;\succeq\;0\quad\text{(Hermitian PSD).}\;} \tag{$\star\!\star\!\star\!\star$}
$$

## 3. Proof

### Step 1. Bochner–Herglotz for $M\,dx$.

$f,g\in[0,1]$ on $[-1,1]$ implies $M(x)=\int f(t)g(x+t)\,dt\ge 0$ for
all $x\in[-2,2]$ (and $M\equiv 0$ outside). So $M\,dx$ is a finite
non-negative Borel measure on $[-2,2]$ (with total mass
$\int_{-2}^{2}M=1$). Bochner's theorem on the interval (equivalently
Herglotz on the dual side, see e.g.\ Akhiezer, *The Classical Moment
Problem*, Ch.\ V §1) gives that the sequence
$\bigl(\hat M(k)\bigr)_{k\in\mathbb Z}$ defined by
$\hat M(k)=\tfrac14\int_{-2}^{2}e^{-i\pi kx/2}M(x)\,dx$ is *positive
definite*, i.e.\ for every finite sequence $u_{0},u_{1},\dots,u_{R}\in\mathbb C$,
$$
\sum_{k,l=0}^{R}u_{k}\overline{u_{l}}\,\hat M(k-l)\;\ge\;0.
$$
Equivalently, the $(R+1)\!\times\!(R+1)$ Hermitian Toeplitz matrix
$\bigl[\hat M(k-l)\bigr]_{0\le k,l\le R}$ is positive semi-definite.

*Proof of the positivity statement* (sketch, included for
self-containedness): the trigonometric polynomial
$p(x):=\bigl|\sum_{k=0}^{R}u_{k}e^{i\pi kx/2}\bigr|^{2}\ge 0$ has the
expansion
$\sum_{k,l}u_{k}\overline{u_{l}}\,e^{i\pi(k-l)x/2}$, so
$\int_{-2}^{2}p(x)M(x)\,dx
=4\sum_{k,l}u_{k}\overline{u_{l}}\,\hat M(k-l)$.
The LHS is $\ge 0$ because $p\ge 0$ and $M\ge 0$. ∎

### Step 2. The program point at admissible $f$ is feasible.

Construct the program point exactly as in v1–v3, with the additional
choice
$$
\tilde A_{k}^{\text{(adm)}}\;:=\;A_{k}\;=\;\frac{4\sin(k\pi/2)}{k\pi}\,a_{k}-2(a_{k}^{2}+b_{k}^{2})
\quad(k=1,\dots,R).
$$

* (A1): by Lemma 5 of White, $A_{k}\in
  [F_{\cos}^{\text{LB}}[k],F_{\cos}^{\text{UB}}[k]]$ (the $\alpha^{\pm}$
  choices satisfy (3.9) since
  $|\!\tfrac{d}{dx}\cos(\pi kx/2)|\le\pi k/2$).
* (A2): equality, by Lemma 3 (3.5).
* $B_{j}^{\text{(adm)}}$: equals the true sine Fourier coefficient
  $-\tfrac{4\sin(j\pi/2)}{j\pi}b_{j}^{\text{true}}$ by Lemma 3 (3.6),
  and the program's $b_{j}^{\text{prog}}$ may be set to
  $b_{j}^{\text{true}}$ via $\delta_{j}=$ full Lemma 4 tail (which
  Lemma 4 itself certifies is in the program's box).

With those choices, the Toeplitz expression $T_{k,l}$ in the program
becomes the *true* $\hat M(k-l)$, so by Step 1
$T_{R+1}\succeq 0$. ✓

### Step 3. Convexity of all new constraints.

* (A1) — two affine inequalities in $(\tilde A_{k},w,v)$. Convex.
* (A2) for **odd** $k$ — $\tilde A_{k}$ below a concave function of
  $(a_{k},b_{k})$: hypograph of concave = convex set in
  $(\tilde A_{k},a_{k},b_{k})$.
* (A2) for **even** $k=2m$ — $\tilde A_{2m}\le-\tfrac12(c_{m}^{2}+d_{m}^{2})$
  is $\tilde A_{2m}\le$ negative-norm-squared, again hypograph of a
  concave function. Convex.
* ($\star\!\star\!\star\!\star$) — the entries of $T_{R+1}$ are
  *affine* functions of the program variables
  $\tilde A_{k},b_{k}^{\text{prog}}$ (the constant $\tfrac14$ on the
  diagonal, $\pm\tfrac12\tilde A_{|k-l|}$ on the real off-diagonals,
  $\pm\tfrac12 B_{|k-l|}$ on the imaginary off-diagonals, with $B$
  itself an affine function of $b^{\text{prog}}$). So $T_{R+1}$ is an
  *affine* matrix-valued function, and
  $T_{R+1}\succeq 0$ is a Linear Matrix Inequality — convex.

Combining all the above with the convex constraints of v3, the
augmented v4 program is convex. ✓

### Step 4. Implementation as a real SDP.

The Hermitian $T_{R+1}\succeq 0$ is equivalent to the real symmetric
constraint
$$
S\;:=\;\begin{pmatrix}\Re T_{R+1} & -\Im T_{R+1}\\ \Im T_{R+1} & \Re T_{R+1}\end{pmatrix}\;\succeq\;0,
$$
a $(2R+2)\!\times\!(2R+2)$ real PSD constraint with affine entries —
a standard LMI handled natively by MOSEK or SCS. With the test
configuration $R=10$ the LMI has size $22\times 22$, modest.

## 4. Why this is genuinely new

Energy gives only the $\ell^{2}$ Bessel-type bound
$\sum_{k\ge1}|\hat M(k)|^{2}\le \tfrac{\Omega-1/4}{2}$. Toeplitz PSD
*also* encodes:

* the $2\times 2$ principal minors:
  $|\hat M(k)|^{2}\le |\hat M(0)|^{2}=\tfrac{1}{16}$ (these are individually weaker than energy but apply uniformly);
* every higher-order principal minor — cubic, quartic, … inequalities
  in $\{\hat M(k)\}_{k=1}^{R}$ that **the energy bound cannot see at
  all**. For example the $3\times 3$ Carathéodory–Toeplitz inequality
  $$\bigl|\hat M(2)-2\hat M(1)^{2}/\hat M(0)\bigr|\le \hat M(0)-\frac{|\hat M(1)|^{2}}{\hat M(0)},$$
  which couples $\hat M(1)$ and $\hat M(2)$ in a strictly non-trivial
  way.

No previous constraint in v0–v3 captures any such correlation between
two different Fourier coefficients of $M$, because energy bounds each
$|\hat M(k)|^{2}$ separately. v4 is the *first* constraint that uses
the **joint distribution** of $\{\hat M(k)\}_{k}$ forced by
$M\ge0$.

## 5. Strictness

For any admissible $f$ the program is feasible (Step 2), so
$\Omega_{\text{v4}}^{\star}\le \inf_{f}\|M\|_{\infty}$ (validity).
Conversely, the v4 feasible set is a *subset* of the v3 feasible set
(it adds a new PSD LMI on top of all v3 constraints; the lifting in
(A1)+(A2) for even $k$ is *consistent* with v3 because $\tilde A_{2m}$
is a new variable that does not appear in v3 and can always be set to
the v3-implicit $-\tfrac12(c_{m}^{2}+d_{m}^{2})$). Therefore
$\Omega_{\text{v4}}^{\star}\ge \Omega_{\text{v3}}^{\star}$.
The new LMI is binding whenever some principal minor of the
v3-optimal Toeplitz matrix is strictly negative — which is generic
unless the v3 optimiser happens to land on a positive-measure-feasible
configuration (in particular, generic in any regime where v3 fails to
fully recover the sign correlations of the true cosine/sine Fourier
coefficients of $M$).

## 6. Note on the program's $b_{k}^{\text{prog}}$ inside $T$

The Toeplitz expression uses the program's truncated representation
$b_{k}^{\text{prog}}=\delta_{k}+\tfrac{4}{\pi}\sum_{j=1}^{T}\tfrac{j(-1)^{j}\sin(k\pi/2)}{k^{2}-4j^{2}}d_{j}$
for the sine Fourier coefficient. As established in the v2 verdict
(§3, "lifting argument"), every admissible $f$ extends to a program
point at which $b_{k}^{\text{prog}}=b_{k}^{\text{true}}$, by choosing
$\delta_{k}=$ Lemma 4 tail. At that admissible point the Toeplitz is
exactly $[\hat M(k-l)]$ and ($\star\!\star\!\star\!\star$) holds by
Bochner. So the use of $b_{k}^{\text{prog}}$ in $T$ is valid.
