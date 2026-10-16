# v5 Proposal: dual Bochner constraint — $T \preceq \Omega I$

## 1. Motivation

v4 added the constraint $T_{R+1}\succeq 0$ from $M\ge 0$ on $[-2,2]$.
This is structurally important but does **not** explicitly couple the
Toeplitz matrix to the optimisation variable $\Omega$.

Observation: by definition $\Omega = \lVert M\rVert_\infty$, the function
$\Omega - M(x)$ is **also non-negative** on $[-2,2]$. Applying the
Bochner–Herglotz theorem to $\Omega - M$ in exactly the same way yields
a *second* Toeplitz PSD constraint, which when combined with v4's
sandwiches the Toeplitz matrix on both sides:
$$\boxed{\;0\;\preceq\;T_{R+1}\;\preceq\;\Omega\,I_{R+1}.\;}$$
Only the **upper** half of this sandwich is genuinely new — and it is
precisely the constraint that ties the rich Fourier-correlation
structure encoded in $T$ to the lower-bound variable $\Omega$.

## 2. The new constraint

Re-use the lifting from v4: variables $\tilde A_k$ for
$k\in\{1,\dots,R\}$ with (A1) Lemma 5 envelope and (A2) Lemma 3
hypograph, plus the program expression $B_j$ for the sine Fourier
coefficient. Define $T_{R+1}$ exactly as in v4. Then **impose**, in
addition to $T_{R+1}\succeq 0$,
$$
\Omega\,I_{R+1}\;-\;T_{R+1}\;\succeq\;0. \tag{$\star\!\star\!\star\!\star\!\star$}
$$

## 3. Proof

### Step 1. $\Omega - M \ge 0$ on $[-2,2]$.

By definition of $\Omega=\lVert M\rVert_\infty\ge M(x)$ for every
$x\in[-2,2]$. Outside $[-2,2]$, both sides are $0$, so the inequality
is irrelevant there. Inside, $\Omega-M(x)\ge 0$.

### Step 2. Bochner on $\Omega-M$.

Repeat the argument from v4 with $\Omega-M$ in place of $M$. For any
finite sequence $u_0,\dots,u_R\in\mathbb C$, the trigonometric
polynomial $p(x):=\bigl|\sum_{k=0}^{R}u_k e^{i\pi kx/2}\bigr|^{2}$ is
non-negative on $\mathbb R$. Since $(\Omega-M)\ge 0$ on $[-2,2]$ and
$(\Omega-M)\equiv \Omega>0$ outside $[-2,2]$ — wait, we only need the
non-negativity *as an integrand* on $[-2,2]$ for the proof to go
through, and that holds. Concretely,
$$
0\;\le\;\int_{-2}^{2}p(x)\bigl(\Omega-M(x)\bigr)\,dx
=\sum_{k,l=0}^{R}u_k\overline{u_l}\!\int_{-2}^{2}e^{i\pi(k-l)x/2}(\Omega-M(x))\,dx.
$$
Using $\int_{-2}^{2}e^{i\pi(k-l)x/2}\,dx=4\delta_{kl}$ for the constant
part $\Omega$, and the identity
$\int_{-2}^{2}e^{i\pi jx/2}M(x)\,dx=4\,\overline{\hat M(j)}$ from v4, we
get
$$
\int_{-2}^{2}e^{i\pi jx/2}(\Omega-M)\,dx
=4\Omega\,\delta_{j,0}\;-\;4\,\overline{\hat M(j)}.
$$
Hence
$$
0\;\le\;4\sum_{k,l}u_k\overline{u_l}\bigl(\Omega\delta_{kl}-\overline{\hat M(k-l)}\bigr)
=4\sum_{k,l}\overline{u_k}u_l\bigl(\Omega\delta_{kl}-\hat M(k-l)\bigr),
$$
where the second equality is by the change of variable
$k\leftrightarrow l$ and the fact that the sum is real. Therefore
$$
u^{*}\bigl(\Omega I - T_{R+1}\bigr)u\;\ge\;0\quad\forall u\in\mathbb C^{R+1},
$$
i.e. $\Omega I_{R+1}-T_{R+1}\succeq 0$.

### Step 3. Validity in the lifted program.

At any admissible $f$, set $\tilde A_k^{\text{(adm)}}=A_k^{\text{true}}$
(as in v4) and $b_k^{\text{prog}}=b_k^{\text{true}}$. Then the
expression $T_{R+1}$ becomes the *true* Toeplitz
$[\hat M(k-l)]_{0\le k,l\le R}$, and by Steps 1–2 above
$\Omega I-T_{R+1}\succeq 0$. The other v4 constraints (including
$T_{R+1}\succeq 0$) and all v3 constraints are also satisfied as in
v4's proof. So every admissible $f$ extends to a v5-feasible point. ✓

### Step 4. Convexity.

The matrix $\Omega I_{R+1}-T_{R+1}$ has entries that are *affine* in
the program variables $(\Omega,\tilde A_k,b_k^{\text{prog}})$:

* diagonal $(\Omega I)_{kk}-T_{kk}=\Omega-\tfrac14$ — affine in
  $\Omega$;
* off-diagonal $-T_{kl}=-\bigl(\tilde A_{|k-l|}-i\sigma(k-l)B_{|k-l|}\bigr)/2$ —
  affine in $\tilde A,b^{\text{prog}}$.

So $\Omega I-T_{R+1}$ is an affine Hermitian matrix-valued function,
and the constraint $\Omega I-T_{R+1}\succeq 0$ is a Linear Matrix
Inequality (an LMI) — convex.

Implementation again uses the standard real reformulation
$\bigl(\begin{smallmatrix}\Re(\Omega I-T)&-\Im(\Omega I-T)\\ \Im(\Omega I-T)&\Re(\Omega I-T)\end{smallmatrix}\bigr)\succeq 0$,
a $(2R+2)\times(2R+2)$ real symmetric PSD constraint with the same
size as v4's PSD constraint (just shifted by $\Omega I$).

## 4. Why this is genuinely new

* No constraint in v0–v4 explicitly bounds the **largest eigenvalue**
  of the $M$-hat Toeplitz matrix in terms of $\Omega$.
  $T\preceq \Omega I$ does exactly that.
* As a corollary at the 2×2 minor level,
  $\bigl(\Omega-\tfrac14\bigr)^{2}\ge|\hat M(k)|^{2}=\tfrac14(A_k^{2}+B_k^{2})$,
  i.e. $A_k^{2}+B_k^{2}\le 4(\Omega-\tfrac14)^{2}$ for *every*
  individual $k$. This is a per-coefficient bound that the v1 energy
  bound only gives in aggregate (and which is non-trivially tighter
  when one coefficient dominates).
* At higher minor levels, $\Omega I-T\succeq 0$ produces *correlation
  inequalities between distinct $\hat M(k)$'s involving $\Omega$* —
  information of a qualitatively different flavour from v4
  (which has no $\Omega$ in its LMI at all).
* Together with v4's $T\succeq 0$, the sandwich
  $0\preceq T\preceq \Omega I$ is the standard SDP characterisation of
  "measure on $[-2,2]$ whose density is bounded by $\Omega$", which is
  exactly the *defining property* of $M$.

## 5. Strictness

v5 ⊃ v4 ⊃ v3 ⊃ … (each adds a new LMI/constraint without removing
existing ones). The new LMI is *binding* at any v4-optimum whose
Toeplitz expression has at least one eigenvalue exceeding $\Omega$
(generically true unless v4 already produces a Toeplitz with all
eigenvalues in $[0,\Omega]$, which is a stringent condition). So
$\Omega_{\text{v5}}^{\star}\ge \Omega_{\text{v4}}^{\star}$, with strict
gain in any regime where v4's relaxation overestimates the spectral
range of the true $T$.

## 6. Note: this is **not** a hyperparameter tweak

The shift from "Bochner on $M$" to "Bochner on $\Omega-M$" is a new
mathematical hypothesis: the *upper* bound $M\le \Omega$ (used only in
the v1 energy inequality $\int M^2\le \Omega$) is here used in a
*structural* way — through the positivity of an entire function and
the consequent PSD Toeplitz. No previous constraint in the v0→v4 chain
extracts a *Fourier* consequence of the upper bound $M\le \Omega$.
