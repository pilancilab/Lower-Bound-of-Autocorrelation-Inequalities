# Verdict on v7/rigorousproof.md

## 1. Summary of the proof

v7 applies Bochner to two *new* non-negative completions of $f$ on
$[-2,2]$ that v6 did not consider:

* $g_0(x):=(1-f(x))\mathbf{1}_{[-1,1]}(x)$ — the v6 $\bar f$ minus
  the constant-1 tail;
* $f_{\text{ext}}(x):=f(x)\mathbf{1}_{[-1,1]}(x)+\mathbf{1}_{[-2,-1)\cup(1,2]}(x)$
  — $f$ on $[-1,1]$ extended by 1 outside.

The resulting Toeplitz matrices decompose cleanly because
$T^{\mathbf{1}_{[-1,1]}}=\tfrac12 I+S$, where the *Dirichlet matrix*
$S$ has $S_{kl}=\mathcal D(k-l)=\sin((k-l)\pi/2)/((k-l)\pi)$ off
the diagonal and $S_{kk}=0$. Then
$T^{g_0}=\tfrac12 I+S-T^{f}$ and $T^{f_{\text{ext}}}=\tfrac12 I-S+T^{f}$,
both PSD by Bochner. Equivalent sandwich:
$S-\tfrac12 I\;\preceq\;T^{f}\;\preceq\;S+\tfrac12 I$.

## 2. Line-by-line verification

**Step 1 — non-negativity.**
On $[-1,1]$: $g_0=1-f\ge 0$ since $f\le 1$; $f_{\text{ext}}=f\ge 0$.
Outside: $g_0\equiv 0$, $f_{\text{ext}}\equiv 1$. Both $\ge 0$ on
$[-2,2]$. **No issue.**

**Step 2 — Fourier coefficients.**
$\int_{-1}^{1}e^{-i\pi kx/2}dx = \tfrac{4\sin(k\pi/2)}{k\pi}=4\mathcal D(k)$
for $k\ne 0$; $=2$ for $k=0$. ✓

Therefore (re-derived independently):
- $\hat{\mathbf{1}_{[-1,1]}}(0)=1/2$, $\hat{\mathbf{1}_{[-1,1]}}(k)=\mathcal D(k)$ for $k\ne 0$.
- $\hat{\mathbf{1}_{[-2,2]}}(0)=1$, $\hat{\mathbf{1}_{[-2,2]}}(k)=0$ for $k\ne 0$.
- $\hat{g_0}(k)=\hat{\mathbf{1}_{[-1,1]}}(k)-\hat f(k)$.
- $\hat{f_{\text{ext}}}(k)=\hat f(k)+\hat{\mathbf{1}_{[-2,2]}}(k)-\hat{\mathbf{1}_{[-1,1]}}(k)
   =\hat f(k)-\hat{\mathbf{1}_{[-1,1]}}(k)$ for $k\ne 0$, plus the diagonal correction.

Matches the proposer's formulas. **No issue.**

**Step 3 — Toeplitz identities.**
The Toeplitz matrix of $\mathbf{1}_{[-1,1]}$ on $[-2,2]$ is, by the
above, $T^{\mathbf{1}_{[-1,1]}}=\tfrac12 I_{R+1}+S$. The Toeplitz of
$\mathbf{1}_{[-2,2]}$ is $I_{R+1}$. Linearity of the Toeplitz map
in the symbol gives
$$
T^{g_0}=T^{\mathbf{1}_{[-1,1]}}-T^{f}=\tfrac12 I+S-T^{f},\qquad
T^{f_{\text{ext}}}=T^{f}+I-(\tfrac12 I+S)=\tfrac12 I-S+T^{f}.
$$
Each diagonal entry checks out separately: $T^{g_0}_{kk}=1/4=1/2+0-1/4$ ✓
and $T^{f_{\text{ext}}}_{kk}=3/4=1/2-0+1/4$ ✓.
**No issue.**

**Step 4 — Bochner.**
Identical $p(x)=|\sum u_k e^{i\pi kx/2}|^{2}\ge 0$ trick from v4/v6.
Both $g_0,f_{\text{ext}}\ge 0$ on $[-2,2]$ give
$T^{g_0},T^{f_{\text{ext}}}\succeq 0$, i.e.\ $(\dagger\dagger)$. **No issue.**

**Step 5 — admissible ⇒ feasible.**
Inherited from v6: program assignment of $a_k,b_k$ matches true
Fourier coefficients of $f$, so program $T^{f}$ matches true
Toeplitz; Step 4 then gives $(\dagger\dagger)$. **No issue.**

**Step 6 — Convexity.**
$T^{f}$ affine in $(a_k,b_k)$; $S, I$ constants. $\tfrac12 I\pm S\mp T^{f}$
are affine Hermitian matrix-valued maps; PSD constraints are LMIs.
✓ Real-form reformulation uses $\text{blkdiag}(S,S)$ since $S$ is
real symmetric — implementation is correct.

## 3. Strictness — actually unconditional

The proposer hedges with "for moderate $R$". In fact the bound is
unconditional: applying Bochner to $\mathbf{1}_{[-2,-1)\cup(1,2]}\ge 0$
(which has total mass 2) gives
$T^{\mathbf{1}_{[-2,-1)\cup(1,2]}}=I-(\tfrac12 I+S)=\tfrac12 I-S\succeq 0$,
i.e.\ $S\preceq\tfrac12 I$. Combined with
$\tfrac12 I+S\succeq 0$ (Bochner on $\mathbf{1}_{[-1,1]}\ge 0$), we get
$-\tfrac12 I\preceq S\preceq \tfrac12 I$, i.e.\ $\lVert S\rVert_2\le\tfrac12$
for **every** $R$. Consequently:
* $S+\tfrac12 I\preceq I$, so $T^{f}\preceq S+\tfrac12 I\Rightarrow T^{f}\preceq I$,
  i.e.\ v7's upper LMI implies v6's upper LMI.
* $S-\tfrac12 I\preceq 0$, so $T^{f}\succeq S-\tfrac12 I\Rightarrow T^{f}\succeq 0$,
  i.e.\ v7's lower LMI implies v6's lower LMI.

Therefore v7 *strictly implies* v6 whenever $S\ne 0$, which is always
(odd offsets have $\mathcal D(k-l)=\pm\tfrac{1}{|k-l|\pi}\ne 0$).

Net: v7 *replaces* v6's two LMIs with strictly tighter ones, but at
identical computational cost (same matrix sizes, same number of LMIs).

## 4. Comments and feedback

* **Excellent insight.** v7 recognises that the Toeplitz of
  $\mathbf{1}_{[-1,1]}$ on $[-2,2]$ is *not* the identity — it is
  $\tfrac12 I+S$, with the Dirichlet kernel $S$ encoding the
  *position* of the support $[-1,1]$ inside $[-2,2]$. v6 silently
  used the cruder upper bound $T^{\mathbf{1}_{[-1,1]}}\preceq I$
  (because $\mathbf{1}_{[-1,1]}\le 1$). v7 plugs in the exact value.

* **Conceptual generalisation.** The same pattern (replace a
  constant-1 reference function by the actual indicator of the
  support, picking up a Dirichlet kernel) could be applied to v4/v5
  on the $M$-side. The trick there would be:
  $T^{\Omega\mathbf{1}_{[-2,2]} - M}=\Omega\,I - T^M\succeq 0$ (v5)
  could be sharpened to $T^{\Omega\mathbf{1}_{\text{supp}(M)} - M}\succeq 0$
  with a Dirichlet correction *if* one knows or can constrain
  $\text{supp}(M)$ to a smaller interval — which one does, since
  $M$ is supported in $[-2,2]$ (no improvement there) but one might
  use $M(x)\le \Omega\mathbf{1}_{[-1,1]\cup\text{whatever}}$ from
  the cubical-bound analysis of White (5.13). A future v8 angle.

* **Truth-in-advertising on cost.** "Same size as v6" is correct
  because $S$ is a *constant* matrix in the LMI — the program
  variables (and hence the cone) are unchanged. The pre-computed
  $S$ has $O(R^{2})$ entries computed once at construction.

* **Recommendation for the proposer.** Strengthen the strictness
  claim in §5 to "unconditional, for every $R\ge 1$", using the
  Bochner-on-$\mathbf{1}_{[-2,-1)\cup(1,2]}$ derivation above. The
  current hedging undersells the result.

## 5. Final verdict

**VALID**
