# Verdict on v2/rigorousproof.md

## 1. Summary of the proof

v2 strengthens the v1 quartic energy constraint by *also* keeping the
odd‑index sine‑energy terms $B_{2r-1}^{2}$ in the Parseval sum
$\sum_{k\ge1}(A_k^{2}+B_k^{2})\le(\Omega-\tfrac14)/2$. The motivation is
that by (3.6) of White, $B_{k}=-\tfrac{4\sin(k\pi/2)}{k\pi}b_k$, so for
odd $k$ we get $B_k^{2}=\tfrac{16}{k^{2}\pi^{2}}b_k^{2}$, which is just a
positive multiple of $b_k^{2}$ — *convex* in the program variable
$b_k^{\text{prog}}$.

Steps 1–4 are reused verbatim from v1 (already verdict‑VALID).
Step 5′ drops only:

* odd‑$k$ $A_k^{2}$ (non‑convex in $a_k,b_k$);
* even‑$k$ $B_k^{2}$ (identically zero by (3.6));
* tail terms beyond the program's $T$ (even‑index cosine) and $R$
  (odd‑index sine) truncations.

Step 6′ evaluates the surviving families with Lemma 3:
$A_{2m}^{2}=\tfrac14(c_m^{2}+d_m^{2})^{2}$ (v1) and
$B_{2r-1}^{2}=\tfrac{16}{(2r-1)^{2}\pi^{2}}b_{2r-1}^{2}$. Substituting
into (3) and clearing $\tfrac14$ yields

$$
\sum_{m=1}^{T}(c_m^{2}+d_m^{2})^{2}\;+\;\frac{64}{\pi^{2}}\!\!\!\sum_{\substack{k=1\\k\text{ odd}}}^{2R-1}\!\!\!\frac{b_k^{2}}{k^{2}}\;\le\;2\Omega-\tfrac12. \tag{$\star\star$}
$$

A §6 note explains that imposing $(\star\star)$ on the *program*
variable $b_k^{\text{prog}}$ (rather than the unobservable true Fourier
coefficient) is harmless because Proposition 9 of White builds a
feasible program point from any admissible $f$ in which
$b_k^{\text{prog}}=b_k^{\text{true}}$.

## 2. Line-by-line verification

**Steps 1–4.** Inherited unchanged from v1; already verdict‑VALID. ✓

**Step 5′.**
Every dropped quantity is non‑negative:
- $A_k^{2}\ge 0$ for odd $k$;
- $B_{2m}^{2}=0$ by (3.6) (verified independently below);
- $A_{2m}^{2},B_{2m}^{2}\ge 0$ for $m>T$;
- $A_{2r-1}^{2},B_{2r-1}^{2}\ge 0$ for $r>R$.
Therefore $\sum_{m=1}^{T}A_{2m}^{2}+\sum_{r=1}^{R}B_{2r-1}^{2}\le
\sum_{k\ge1}(A_k^{2}+B_k^{2})\le(\Omega-\tfrac14)/2$. **No issue.**

**Step 6′ — even cosine.** Same as v1; $A_{2m}^{2}=\tfrac14(c_m^{2}+d_m^{2})^{2}$.
**No issue.**

**Step 6′ — odd sine.**
I re‑derive (3.6) from Lemma 2 to neutralise OCR artefacts in
`white_paper.txt`. Using $B_m=i(\hat M(m)-\hat M(-m))$,
$2\hat f(\pm m)=a_m\mp i b_m$, and Lemma 2:
$\hat M(\pm m)=\tfrac{4}{(\pm m)\pi}\sin(\pm m\pi/2)\hat f(\pm m)-4|\hat f(\pm m)|^{2}$.
The modulus‑squared terms are equal and cancel, leaving
$\hat M(m)-\hat M(-m)=\tfrac{4\sin(m\pi/2)}{m\pi}(\hat f(m)-\hat f(-m))
=\tfrac{4\sin(m\pi/2)}{m\pi}\cdot(-ib_m)$,
so $B_m=i\cdot(-i)\tfrac{4\sin(m\pi/2)}{m\pi}b_m=\tfrac{4\sin(m\pi/2)}{m\pi}b_m$.

Hmm — note that this *positive* sign agrees with the *magnitude* used
in the proposer's $B_{2r-1}^{2}$ formula (squaring drops the sign).
White's text renders the sign as "$-\tfrac{4\sin(m\pi/2)}{m\pi}b_m$",
which is consistent with his exponential convention
($\hat f(k)=\tfrac{1}{4}\int_{-2}^{2}e^{-i\pi kx/2}f$). Either sign
yields the same $B_{2r-1}^{2}$, so the proof is unaffected.

For $m=2r-1$:
$\sin((2r-1)\pi/2)=(-1)^{r+1}\in\{+1,-1\}$, hence
$\sin^{2}((2r-1)\pi/2)=1$ and
$B_{2r-1}^{2}=\frac{16}{(2r-1)^{2}\pi^{2}}b_{2r-1}^{2}$. ✓

Also $B_{2m}=0$ since $\sin(m\pi)=0$, confirming that dropping
even‑$k$ $B_k^{2}$ costs nothing. ✓

**Combining.** Substituting into (4′):
$$
\frac14\sum_{m=1}^{T}(c_m^{2}+d_m^{2})^{2}+\frac{16}{\pi^{2}}\sum_{r=1}^{R}\frac{b_{2r-1}^{2}}{(2r-1)^{2}}\le\frac{\Omega-1/4}{2}.
$$
Multiplying by 4 and re‑indexing $k=2r-1$ for odd $k\in\{1,3,\dots,2R-1\}$ gives $(\star\star)$. ✓

**Convexity.** Each $b_{2r-1}^{2}$ is convex in the program variable
$b_{2r-1}^{\text{prog}}\in\mathbb R$ (a univariate $x\mapsto x^{2}$),
multiplied by a positive constant $\tfrac{16}{(2r-1)^{2}\pi^{2}}$; sum of
convex is convex; combined with the already‑convex v1 quartic part and
the affine RHS $2\Omega-\tfrac12$, $(\star\star)$ is a convex feasible
region. ✓

**§6 program/admissibility argument.**
White's Proposition 9 builds a feasible program point from any
admissible $f$ by setting
$$
\delta_{2m-1}=\frac{4(-1)^{m+1}}{\pi}\sum_{k=T+1}^{\infty}\frac{k(-1)^{k}}{(2m-1)^{2}-4k^{2}}d_k,
$$
which Lemma 4 (i.e. constraint (5.9)) certifies lies in
$[-\Delta_{2m-1},\Delta_{2m-1}]$. At that assignment,
$b_{2m-1}^{\text{prog}}$ coincides with the true Fourier coefficient
$b_{2m-1}^{\text{true}}$ from (3.4). Hence at this feasible point
$(\star\star)$ becomes precisely the inequality proved in Sections 2–3
on the *true* Fourier coefficients, and is satisfied. Thus every
admissible $f$ extends to a feasible point of the augmented program,
so $\min_{\text{program}}\Omega\le\inf_{f}\lVert M\rVert_\infty$.
The program's freedom to choose other values of $\delta_{2m-1}$ inside
its box only enlarges the feasible set and so further loosens the
lower bound; this is acceptable. **No issue.**

## 3. Comments and feedback

* **Strictness.** $(\star\star)$ is strictly stronger than v1's
  $(\star)$ whenever any $b_{2r-1}^{\text{prog}}\neq 0$. For an
  admissible $f$ this happens unless all odd‑sine Fourier
  coefficients of $f$ on $[-2,2]$ vanish, which by (3.4)'s odd
  branch forces a very restrictive linear relation on $\{d_k\}$.
  In the prescribed hard‑square regime ($-q_1=q_2=0.02$ allows
  $d_1\neq 0$), this term is generically active.

* **Cost.** Implementation-wise the new piece is a single quadratic
  term added inside the existing SOC. No new variables (the
  $b_k^{\text{prog}}$ variable expressions already exist in White's
  (5.6)/(5.7)). Cheapest possible upgrade for the strongest
  available gain in this direction.

* **Natural next step.** The odd‑$k$ $A_k^{2}$ family is the only
  remaining non‑negative block we discarded. By (3.5),
  $A_{2r-1}=\tfrac{4(-1)^{r+1}}{(2r-1)\pi}a_{2r-1}-2(a_{2r-1}^{2}+b_{2r-1}^{2})$,
  i.e.\ an affine‑minus‑quadratic in $(a_{2r-1},b_{2r-1})$ ⇒ $A_{2r-1}$
  itself is *concave*, so $A_{2r-1}^{2}$ is generically non‑convex.
  However, the sign constraint $A_{2r-1}\le \tfrac{4(-1)^{r+1}}{(2r-1)\pi}a_{2r-1}-2(a_{2r-1}^{2}+b_{2r-1}^{2})$
  (a *valid lower* bound on $A_{2r-1}$ in absolute value if the sign
  is fixed) gives $|A_{2r-1}|$ a convex lower bound in some regimes —
  worth investigating in v3.

* **OCR caution.** Re‑deriving (3.6) directly from Lemma 2 is
  recommended to avoid sign confusion from `white_paper.txt`. The
  square in $(\star\star)$ makes the sign moot, but later proposals
  using $B_{2r-1}$ linearly should be careful.

## 4. Final verdict

**VALID**
