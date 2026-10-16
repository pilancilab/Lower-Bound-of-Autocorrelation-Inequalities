# Verdict on v4/rigorousproof.md

## 1. Summary of the proof

v4 leaves the energy programme of v1–v3 intact and adds a single LMI:
the $(R+1)\times(R+1)$ Hermitian Toeplitz matrix
$$
T_{R+1}\;:=\;[\hat M(k-l)]_{0\le k,l\le R}
\;=\;
\begin{cases}
1/4 & k=l\\
(\tilde A_{|k-l|}-i\sigma(k-l)B_{|k-l|})/2 & k\neq l
\end{cases}
$$
must be PSD. The constraint is grounded in Bochner–Herglotz:
$M\ge 0$ on $[-2,2]$ ⇒ $(\hat M(k))_{k\in\mathbb Z}$ is a
positive‑definite sequence. To make this work in the convex program,
v4 extends the v3 lifting of $\tilde A_k$ to **all** $k\in\{1,\dots,R\}$
(odd *and* even), each governed by the two convex one-sided
inequalities (A1) (Lemma 5 envelope) and (A2) (Lemma 3 (3.5)
hypograph). $B_{|k-l|}$ is the program's affine expression for the
sine Fourier coefficient (Lemma 3 (3.6) combined with the
Lemma 4 truncation of $b_k^{\text{prog}}$, as already validated in
v2's verdict).

The PSD constraint is the first proposal in the v0→v4 chain that
captures *correlations between distinct Fourier coefficients of $M$*
— every energy bound so far only constrained each $|\hat M(k)|^{2}$
separately.

## 2. Line-by-line verification

**Step 1 — Bochner–Herglotz on a finite interval.**

The argument is the textbook one. Let
$p(x):=|\sum_{k=0}^R u_k e^{i\pi kx/2}|^2\ge 0$. Then
$\int_{-2}^{2} p(x)M(x)\,dx\ge 0$ because both factors are
non-negative on $[-2,2]$. Expanding,
$$
\int_{-2}^{2}p(x)M(x)\,dx
=\sum_{k,l}u_k\overline{u_l}\!\int_{-2}^{2}\!e^{i\pi(k-l)x/2}M(x)\,dx
=4\sum_{k,l}u_k\overline{u_l}\,\overline{\hat M(k-l)},
$$
using White's convention $\hat M(j)=\tfrac14\!\int_{-2}^{2}e^{-i\pi jx/2}M$.
The sum is real (LHS is real and the change of variables
$k\leftrightarrow l$ conjugates each summand), so equivalently
$\sum_{k,l}\overline{u_k}u_l\hat M(k-l)\ge0$, i.e.
$u^{*}T_{R+1}u\ge0$ for all $u\in\mathbb C^{R+1}$. Hence
$T_{R+1}\succeq 0$. **No issue.** The remark about needing
"Bochner on the interval" is appropriate; no global periodic
extension of $M$ is required because we integrate only on $[-2,2]$,
where $M\ge 0$ holds.

**Toeplitz entries match $\hat M$.**
Diagonal: $\hat M(0)=\tfrac14\!\int M=\tfrac14$. ✓
For $j\ne 0$:
$A_{|j|}=\hat M(|j|)+\hat M(-|j|)$ and $B_{|j|}=i(\hat M(|j|)-\hat M(-|j|))$
give $\hat M(|j|)=(A_{|j|}-iB_{|j|})/2$ and
$\hat M(-|j|)=(A_{|j|}+iB_{|j|})/2$. Therefore
$\hat M(k-l)=(A_{|k-l|}-i\sigma(k-l)B_{|k-l|})/2$,
matching the proposer's $T_{k,l}$ when $\tilde A=A^{\text{true}}$. ✓
Hermitian: $T_{lk}=\overline{T_{kl}}$ by direct check. ✓

**$B_j$ formula.**
White's (3.6) gives $B_m=-\tfrac{4\sin(m\pi/2)}{m\pi}b_m$. For
$m=2r-1$, $\sin((2r-1)\pi/2)=(-1)^{r+1}$, so
$B_{2r-1}=-\tfrac{4(-1)^{r+1}}{(2r-1)\pi}b_{2r-1}^{\text{prog}}$;
for $m=2r$, $\sin(r\pi)=0$ ⇒ $B_{2r}=0$. ✓
The use of $b^{\text{prog}}$ vs.\ $b^{\text{true}}$ is reconciled
exactly as in v2's verdict: Lemma 4 + Proposition 9 of White give
$b_k^{\text{prog}}=b_k^{\text{true}}$ at the admissible-extension
point.

**Step 2 — admissible $f$ ⇒ program-feasibility.**

At an admissible $f$, set $\tilde A_k^{\text{(adm)}}:=A_k^{\text{true}}$
for all $k=1,\dots,R$.

* **(A1).** Lemma 5 (verified for v3 via the $\alpha^\pm$ envelope
  satisfying (3.9) by the derivative bound
  $|\!\tfrac{d}{dx}\cos(\pi kx/2)|\le\pi k/2$). ✓ The argument is the
  same for *every* $k$, not only odd.
* **(A2).** Lemma 3 (3.5) is an equality:
  - odd $k=2r-1$:
    $A_k=\tfrac{4(-1)^{r+1}}{k\pi}a_k-2(a_k^2+b_k^2)$;
  - even $k=2m$:
    $A_{2m}=-2(a_{2m}^{2}+b_{2m}^{2})=-\tfrac12(c_m^{2}+d_m^{2})$
    via $a_{2m}=\tfrac12c_m,\ b_{2m}=\tfrac12d_m$
    (Lemma 3 (3.3)/(3.4) even branch).
  Either way (A2) holds with equality. ✓
* **PSD.** With $\tilde A=A^{\text{true}}$ and
  $b^{\text{prog}}=b^{\text{true}}$, $T_{R+1}$ becomes
  $[\hat M(k-l)]$, which is PSD by Step 1. ✓

**No issue.**

**Step 3 — Convexity.**

* (A1): two affine inequalities in
  $(\tilde A_k,\{w_j,v_j\})$. Convex. ✓
* (A2), odd $k$: $\tilde A_k$ below an affine-minus-norm-squared
  function of $(a_k,b_k)$ — concave. Hypograph of concave function
  is a convex set in $(\tilde A_k,a_k,b_k)$. ✓
* (A2), even $k=2m$: $\tilde A_{2m}\le-\tfrac12(c_m^{2}+d_m^{2})$
  is the hypograph of a concave function of $(c_m,d_m)$. Convex. ✓
* $T_{R+1}\succeq 0$ where entries are *affine* functions of the
  program variables $\tilde A_k$ and $b_k^{\text{prog}}$ (the
  scalar $\tfrac12$ multiplier and $\pm i\sigma$ factor are
  constants, $\tilde A_{|k-l|}$ enters linearly, $B_{|k-l|}$ is
  $-\tfrac{4(-1)^{r+1}}{j\pi}\cdot b_j^{\text{prog}}$ which is
  affine in $b_j^{\text{prog}}$ which is itself affine in the
  program's $\delta_j$ and $d_j$). PSD constraint on an affine
  Hermitian matrix-valued map is an LMI. Convex. ✓

**Implementation.** The real-block reformulation
$\Bigl(\begin{smallmatrix}\Re T & -\Im T\\ \Im T & \Re T\end{smallmatrix}\Bigr)\succeq 0$
is the standard Hermitian-to-real-symmetric trick. For $R=10$ the
LMI is $22\times 22$ with affine entries — well within MOSEK/SCS
limits. ✓

## 3. Strictness vs. v3

v4 inherits all of v3 (variables, constraints, energy inequality)
and adds:
1. fresh variables $\tilde A_{2m}$ for $m=1,\dots,\lfloor R/2\rfloor$
   with (A1) and the even-$k$ (A2);
2. the PSD LMI.

(1) is consistent with v3 because each $\tilde A_{2m}$ is brand new
and may be set to the v3-implicit value
$-\tfrac12(c_m^{2}+d_m^{2})$ to recover an admissible v3 point.
(2) is a strict additional constraint. Hence v4-feasible ⊂
v3-feasible and $\Omega_{\text{v4}}^{\star}\ge\Omega_{\text{v3}}^{\star}$.
The PSD constraint is **binding** whenever the v3 optimiser produces
a Toeplitz expression whose smallest eigenvalue is < 0; this is
generic in regimes where v3's $\tilde A_{2r-1}$ choices do not align
with a true positive‑measure-feasible configuration. ✓

## 4. Comments and feedback

* **Effectiveness — high.** This is the first proposal that
  taps a *qualitatively different* property of $M$ ($M\ge 0$,
  via Bochner) rather than yet another shadow of
  $\int M^{2}\le \Omega$. Among the standard moment/dual-cone
  techniques for proving lower bounds on $\|M\|_\infty$, the
  Toeplitz PSD hierarchy is the canonical strengthening; v4
  imports it cleanly using the v3 lifting infrastructure.

* **Hierarchy.** The PSD constraint with size $R+1$ is the
  *first* level of the Carathéodory–Toeplitz / Schur–Cohn moment
  hierarchy. Higher levels (or interlacing with additional
  Fourier blocks) are natural future directions. Note in
  particular that diagonal entries $\hat M(0)=\tfrac14$ are
  *constants*: this hard‑wires the average constraint and gives
  the LMI a unit‑normalised scale.

* **Computational note.** For the validation regime ($R\approx 10$
  by default), the LMI is small — but mixing LMIs with SOCPs has
  historically been a tricky DCP construction in CVXPY. The
  proposer should verify that MOSEK/SCS accepts the augmented
  problem; if SCS struggles, the proof remains valid but
  numerically one may need to drop to a smaller $R$. This is an
  implementation concern, not a proof concern.

* **Tightness of (A2) even-$k$ branch.** $\tilde A_{2m}\le
  -\tfrac12(c_m^{2}+d_m^{2})$ is the *correct* convex side
  (hypograph of concave). The reverse inequality
  $\tilde A_{2m}\ge -\tfrac12(c_m^{2}+d_m^{2})$ is the epigraph
  of a concave function — non-convex — and would be a trap.
  v4 correctly takes only the convex side.

* **Bochner caveat re-checked.** The classical Bochner–Herglotz
  theorem strictly speaking concerns positive-definite functions
  on $\mathbb R$ (Bochner) or sequences of integers (Herglotz). The
  finite-interval statement is implied: it suffices that
  $p(x)\ge 0$ on $[-2,2]$ for the integral inequality. The
  proposer's parenthetical citation of Akhiezer is appropriate
  context but not strictly necessary; the direct argument given
  in Step 1 is self‑contained.

## 5. Final verdict

**VALID**
