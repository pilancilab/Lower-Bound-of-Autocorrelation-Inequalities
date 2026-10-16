# Verdict on v6/rigorousproof.md

## 1. Summary of the proof

v6 imports the Bochner sandwich from v4/v5 to the *input* side. Since
$0\le f\le 1$ on $[-1,1]$ (and $f\equiv 0$ off $[-1,1]$ by convention),
both $f\ge 0$ and the "complemented–and–unit-extended" function
$$
\bar f(x)\;:=\;\begin{cases}1-f(x)&x\in[-1,1]\\ 1&x\in[-2,-1)\cup(1,2]\end{cases}
$$
are non-negative on $[-2,2]$. Applying Bochner to each yields
$$
0\;\preceq\;T_{R+1}^{f}\;\preceq\;I_{R+1},\qquad T_{kl}^{f}:=\hat f(k-l).
$$
Diagonal $T^f_{kk}=\hat f(0)=\tfrac14$ (from $\int_{-1}^{1}f=1$);
off-diagonal $T^f_{kl}=(a_{|k-l|}-i\sigma(k-l)b_{|k-l|})/2$ via
$\hat f(m)=(a_m-ib_m)/2$ (White line 255). The two new LMIs are
imposed on top of all v0–v5 constraints.

## 2. Line-by-line verification

**Step 1 — non-negativity.** $f\ge 0$ on $[-2,2]$ by extension-by-zero
and admissibility. $\bar f\ge 0$ on $[-2,2]$ since $1-f\ge 0$ on
$[-1,1]$ and $\bar f=1$ outside. ✓

**Step 2 — $\hat f(0)$ and $\hat{\bar f}(0)$.**
$\hat f(0)=\tfrac14\!\int_{-2}^{2}f=\tfrac14\!\int_{-1}^{1}f=\tfrac14$.
$\hat{\bar f}(0)=\tfrac14[\int_{-1}^{1}(1-f)+\int_{[-2,-1)\cup(1,2]}1]
=\tfrac14[1+2]=\tfrac34$. ✓

**Step 3 — $\hat{\bar f}(k)$ for $k\ne 0$.**
Independent derivation: $\bar f = \mathbf 1_{[-2,2]} - f$ (with $f$
extended by 0). For integer $k\ne 0$,
$\hat{\mathbf 1_{[-2,2]}}(k)=\tfrac14\!\int_{-2}^{2}e^{-i\pi kx/2}dx=0$
(period-4 exponential over its full period). So
$\hat{\bar f}(k)=0-\hat f(k)=-\hat f(k)$, matching the proposer. ✓

Combining diagonal and off-diagonal:
$T^{\bar f}_{kl}=\hat{\bar f}(k-l)=\delta_{kl}\cdot\tfrac34-(1-\delta_{kl})\hat f(k-l)
=\delta_{kl}(1-\tfrac14)-(1-\delta_{kl})\hat f(k-l)
=I_{R+1}-T^{f}_{R+1}$
(the off-diagonals of $I$ are 0, contributing $-\hat f$; the diagonal
is $1-\tfrac14=\tfrac34$). ✓

**Step 4 — Bochner.** Same trigonometric-polynomial argument as v4:
$p(x)=|\sum u_k e^{i\pi kx/2}|^{2}\ge 0$, both $f\ge 0$ and
$\bar f\ge 0$, so the integrals
$\int_{-2}^{2}p\cdot f$ and $\int_{-2}^{2}p\cdot \bar f$ are $\ge 0$,
which translates respectively to $T^{f}_{R+1}\succeq 0$ and
$T^{\bar f}_{R+1}=I-T^{f}_{R+1}\succeq 0$. ✓

**Step 5 — Convexity.** Entries of $T^{f}_{R+1}$ are *affine* in the
program variables $a_k,b_k$ (the diagonal is the *constant* $1/4$).
$T^{f}$ is Hermitian. Both $T^{f}\succeq 0$ and $I-T^{f}\succeq 0$
are LMIs on an affine Hermitian matrix-valued map — convex. Real
reformulation is the standard Hermitian-to-real-symmetric trick at
size $2(R+1)$. ✓

**Step 6 — admissible ⇒ feasible.** Program assignment for $a_k,b_k$
inherited from v0–v5 (Lemma 4 tail slacks set so
$a_k^{\text{prog}}=a_k^{\text{true}}$ and similarly for $b_k$). At
this assignment $T^f$ in the program equals the true Toeplitz
$[\hat f(k-l)]$, and Bochner (Step 4) gives both LMIs. All other
constraints already verified in v0–v5. ✓

## 3. Independence from existing constraints

**Independence from (5.10) and (5.11).** (5.10) bounds $|c_k|,|d_k|\le 2/\pi$;
(5.11) gives $\sum(c_k^{2}+d_k^{2})\le 1/2$. Neither sees off-diagonal
*correlations* of $T^{f}$: e.g.\ the $2\times 2$ minor of $T^{f}\succeq 0$
on indices $(0,k)$ gives $|\hat f(k)|^{2}\le\tfrac14\cdot\tfrac14=\tfrac{1}{16}$,
i.e.\ $a_k^{2}+b_k^{2}\le\tfrac14$ for *each individual* $k$ — a
per-coordinate bound the energy aggregate cannot enforce. Higher-order
minors give cubic, quartic, … correlations between distinct
$\hat f(k_1),\hat f(k_2),\dots$. The $I-T^{f}\succeq 0$ half is the
$f\le 1$ side and is similarly independent. ✓

**Independence from White's Lemma 3 identities.** Lemma 3 gives exact
linear/quadratic *equalities* tying $(a,b)$ to $(c,d,A,B)$ — no
positivity. v6 is the first constraint imposing the *Bochner*
positivity directly on $(a,b)$. ✓

## 4. Strictness vs.\ v5

v6 inherits all v5 constraints and adds the two LMIs $(\dagger)$.
v6-feasible ⊂ v5-feasible, so
$\Omega_{\text{v6}}^{\star}\ge\Omega_{\text{v5}}^{\star}$, with strict
gain whenever the v5 optimum's induced $T^{f}$ has an eigenvalue
outside $[0,1]$. Since v0–v5 only enforce $L^{p}$ shadows of $f\in[0,1]$,
generically the v5 optimum's $(a,b)$ does not correspond to any actual
non-negative density bounded by 1, so the gap is real. ✓

## 5. Comments and feedback

* **Symmetry with v4/v5.** Conceptually this is the perfect dual of
  the v4+v5 pair. v4+v5 give the SDP characterisation of "$M$ is a
  measure with bounded density in $[0,\Omega]$" truncated to $R+1$
  Fourier modes. v6 gives the SDP characterisation of "$f$ is a
  measure with bounded density in $[0,1]$" truncated to $R+1$
  Fourier modes. After v6 the Toeplitz programme is *symmetric*
  across input and output sides.

* **No extra lifting needed.** Unlike v3/v4 (which had to introduce
  fresh $\tilde A_k$ variables because the Lemma 3 expressions were
  concave-in-program-variables), v6 uses $a_k,b_k$ *directly*: they
  are already linear program variables (or affine expressions of the
  $c_k,d_k$ and $\varepsilon_k,\delta_k$ slacks). The LMIs are
  immediately well-posed.

* **Cheapness.** Each new LMI is $(2R+2)\times(2R+2)$ in real form,
  same size as v4/v5. A solver that already accepts v5 will accept
  v6 with no significant overhead.

* **Possible cross-coupling for future versions.** The next natural
  step is a *joint* moment matrix coupling $\hat f$ and $\hat M$
  (e.g.\ exploiting $\hat M$-via-Lemma-2 as an affine relation between
  the two Toeplitz blocks). This would push beyond the
  "marginal Bochner" content of v4–v6 into genuinely joint structure.

* **Caveat re $\bar f$ extension.** The proposer wisely chose
  $\bar f=1$ outside $[-1,1]$, *not* the periodic extension of
  $1-f$. The periodic extension of $1-f$ from $[-1,1]$ to $[-2,2]$
  would NOT be non-negative outside $[-1,1]$ (it would simply equal
  $1-f$ shifted, which isn't defined the same way). The constant-1
  extension is the *correct* non-negative completion. ✓

## 6. Final verdict

**VALID**
