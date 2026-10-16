# Verdict — Autocorrelation Inequality 6.2, v50

## 1. Summary

v50 introduces a **genuinely new PSD constraint**: the Gram matrix
\(G\) of the real trigonometric basis \(\{1, \cos(2\pi k x), \sin(2\pi k x)\}_{k=1}^{K}\)
in the \(f\)-weighted inner product
\(\langle \phi, \psi\rangle_f = \int f\,\phi\,\psi\,dx\), extended
via auxiliary variables \(a_{K+1..2K}, b_{K+1..2K}\):
\[
G[i, j] \;=\; \int f\,g_i\,g_j\,dx, \quad 0\le i, j\le 2K,
\]
with entries linear in \((a_k, b_k)_{k=0}^{2K}\) via the standard
trig product-to-sum identities. The constraint is
\(G\succeq 0\).

v50 = v45 + (aux \(a, b\) extension from v46) + \(G\succeq 0\).

No numerical results yet.

## 2. Verdict

**VALID** — new and substantive construction, soundly derived.

## 3. Verification

### Gram-matrix entries. ✓

Using the standard trig product-to-sum identities:
\[
\cos(2\pi i x)\cos(2\pi j x) = \tfrac{1}{2}\bigl[\cos(2\pi|i-j|x) + \cos(2\pi(i+j)x)\bigr]
\]
\[
\sin(2\pi i x)\sin(2\pi j x) = \tfrac{1}{2}\bigl[\cos(2\pi|i-j|x) - \cos(2\pi(i+j)x)\bigr]
\]
\[
\cos(2\pi i x)\sin(2\pi j x) = \tfrac{1}{2}\bigl[\sin(2\pi(i+j)x) - \sin(2\pi(i-j)x)\bigr]
\]

Integrating against \(f\) and using \(a_k = \int f\cos(2\pi k x)\),
\(b_k = \int f\sin(2\pi k x)\):

* \(G[0, 0] = 1\), \(G[0, i] = a_i\), \(G[0, K+j] = b_j\). ✓
* \(G[i, j] = \tfrac{1}{2}(a_{|i-j|} + a_{i+j})\) (cos-cos block). ✓
* \(G[K+i, K+j] = \tfrac{1}{2}(a_{|i-j|} - a_{i+j})\) (sin-sin). ✓
* \(G[i, K+j] = \tfrac{1}{2}(b_{i+j} + \text{sgn}(j - i)\,b_{|j-i|})\)
  (cos-sin), with \(\text{sgn}(0) = 0\) giving \(G[i, K+i] = b_{2i}/2\). ✓

The cross-sign factor \(\text{sgn}(j - i)\) correctly accounts for
\(\sin\) being odd. ✓

All entries involve \((a_k, b_k)\) for \(k = 0, 1, \dots, 2K\), so
the auxiliary variables from v46's extension are required.

### Gram PSD. ✓

For any \(v\in\mathbb R^{2K+1}\), with \(\phi(x) = \sum_i v_i g_i(x)\),
\[
v^{T}G v \;=\; \sum_{i, j}v_i v_j \int f g_i g_j\,dx \;=\; \int f(x)\,\phi(x)^{2}\,dx \;\ge\; 0.
\]
So \(G\succeq 0\) for every nonneg \(f\). ✓

### Validity of the relaxation.

At every rank-one truth from admissible \(f\), the Fourier moments
\((a_k, b_k)_{k=0}^{2K}\) satisfy the trig-identity formulas for
\(G\), and \(G\succeq 0\) by the argument above. So the rank-one
truth is feasible in v50, and adding \(G\succeq 0\) to the LP is a
valid relaxation. As a minimization, the LP optimum can only
**increase**: \(\Omega^{(50)}\ge \Omega^{(45)}\). ✓

### Why \(G\) is a genuinely new constraint.

This is the substantive observation:

1. **v14/v45's moment matrix \(M\)** is \((2K+1)\times(2K+1)\) on
   the basis \((1, a_1, \dots, a_K, b_1, \dots, b_K)\), with
   **free SDP entries** \(M[i, j]\) constrained only by
   \(M[0, k] = a_k\), \(M[0, K+k] = b_k\), and \(M\succeq 0\).
   The interior entries \(M[i, j]\) for \(i, j\ge 1\) are SDP
   variables.
2. **v50's Gram matrix \(G\)** has **entries linearly determined**
   by \((a_k, b_k)_{k=0}^{2K}\) via the trig identities. The
   entries are not free; they are affine in the (extended) Fourier
   variables.

At the rank-one truth, \(M[i, j]\) "should" equal \(G[i, j]\) —
both are just \(\int f g_i g_j\,dx\). But in the relaxation, \(M\)
is free while \(G\) is pinned to the linear structure. Adding
\(G\succeq 0\) imposes a **new constraint** on the extended Fourier
sequence, independent of \(M\).

In particular, \(G\) involves \(a_{i+j}\) for \(i + j\) up to
\(2K\) (via the product-to-sum). The Bochner Toeplitz \(T\) only
involves \(a_{|i-j|}\). So \(G\) provides positivity information
at a *different* set of Fourier moments than Bochner.

### Relationship to v47's Cauchy–Schwarz.

v47 gave the diagonal constraints
\(M[k, k]\le (1 + a_{2k})/2 = G[k, k]\)
and \(M[K+k, K+k]\le (1 - a_{2k})/2 = G[K+k, K+k]\). These are
exactly the diagonal entries of \(G\). But v47 constrains \(M\)
(the moment lift), not the Gram matrix.

v50's \(G\succeq 0\) is a PSD constraint on a **different** matrix
(the rank-1-linear Gram, not the SDP-free \(M\)). The two are
complementary, not equivalent. Adding both v47 and v50 may give
further tightening.

## 4. Comments

### Comment 1. This is the right kind of constraint.

After many iterations of vacuous tightening attempts (v15–v28,
v32, v36, v40), v50 is the first *genuinely new* PSD constraint on
the Fourier sequence that goes beyond truncated Bochner. The
\(G\) matrix captures the rank-1 structure of the moment matrix
in the trig basis, which the SDP-free \(M\) does not automatically
impose.

The expected gain depends on whether the SDP's \(M\) at the optimum
is already close to rank-1 (in which case \(G\succeq 0\) is
vacuous) or whether \(M\) exploits the free-entry freedom (in which
case \(G\succeq 0\) tightens). The proposer's small-parameter
sanity run should reveal which.

### Comment 2. The product-to-sum coupling is the key.

\(G[i, j]\) couples low-frequency variables \(a_{|i-j|}\) with
*high-frequency* variables \(a_{i+j}\) (up to \(2K\)). This is the
structural feature that makes \(G\) non-trivial. The low-frequency
entries can be driven to the LP's advantage, but the high-frequency
entries are bounded by cell rectangles, and \(G\succeq 0\) enforces
a PSD coupling between them.

### Comment 3. Expected gain — uncertain but potentially large.

At small parameters where the v45 \(M\) lift is already tight,
\(G\succeq 0\) will give a small gain. At larger \(K\) where the
\(M\) lift has more freedom, \(G\succeq 0\) may give a
substantially larger gain.

### Comment 4. Inheritance is clean.

v50 inherits from v45 + v46 auxiliary variables. No new bugs.

## 5. Feedback to the proposer

For v51:

1. **Run v50 at small-parameter sanity check**: \((8, 8), (16, 16), (32, 32), (50, 32)\).
   Verify \(\Omega^{(50)}\ge \Omega^{(45)}\).
2. **Run v50 at production parameters**: \((1000, 64), (1000, 96), (3000, 96), (1000, 128)\).
   Compare to v44 and v45.
3. **If v50 gives a substantial gain** (\(\ge 10^{-3}\)), the new
   construction is the project's best, and should be combined with
   v47 and v48's \(K\)-extension.
4. **If v50 is vacuous** (\(\le 10^{-5}\)), the SDP's \(M\) lift
   is already effectively rank-1 at the optimum, and the
   framework's ceiling is at \(\sim 1.294\).

---

**Final verdict: VALID** (new and substantive: the Gram matrix
\(G\) of the real trig basis in the \(f\)-weighted inner product,
with entries linear in the extended Fourier variables via
product-to-sum identities, constrained to be PSD. The construction
is correctly derived, and \(\Omega^{(50)}\ge \Omega^{(45)}\)
provably. Numerical results pending. Best bound on \(C_{6.2}\)
remains v44's \(1.2913\) at \((3000, 96)\)).

---

## 6. Independent re-review (skeptical second pass, 2026-04-14) — VERDICT DOWNGRADED

**The earlier verdict is wrong as a claim of tightening when v46's
Hermitian-Toeplitz-at-size-\(2K+1\) is already enforced, as it is
in v50's own code.** Revised verdict:

**INVALID AS A TIGHTENING (tautologically redundant)** given
v46's size-\(2K+1\) Bochner–Toeplitz.

### Why

The real-basis Gram \(G\) in v50 is related to the complex-basis
Hermitian Toeplitz \(T\) at size \(2K+1\) (with frequencies
\(-K, \dots, K\)) by an **invertible congruence**:
\[
G \;=\; P^{*}\,T\,P,
\]
where \(P\in\mathbb C^{(2K+1)\times(2K+1)}\) is the (non-unitary
but invertible) change-of-basis from the complex exponentials to
the real trig basis \(\{1,\cos_k,\sin_k\}\). Explicitly:
* column \(0\) of \(P\) is \(e_{0}\),
* column \(k\) (\(1\le k\le K\)) places \(1/2\) at rows \(\pm k\),
* column \(K+k\) places \(+i/2\) at row \(-k\) and \(-i/2\) at row \(+k\).

Because \(P\) is invertible, PSD is preserved *both ways*:
\(T\succeq 0\Leftrightarrow G\succeq 0\). Hence
\(\{G\succeq 0\}\) and \(\{T^{(2K+1)}\succeq 0\}\) cut the **same**
region of the moment space \((a_0,\dots,a_{2K},b_1,\dots,b_{2K})\).

Centered Toeplitz on \(\{-K,\dots,K\}\) is the same matrix (after
relabelling rows by \(+K\)) as the Toeplitz on \(\{0,\dots,2K\}\)
that v46 / v50's code actually enforces at
`convex_code.py:106–128`. So v46's constraint in v50's code is
exactly \(T\succeq 0\), and \(G\succeq 0\) adds nothing.

### Numerical confirmation

I ran 1000 random draws of \((a_0,\dots,a_4,b_1,\dots,b_4)\) at
\(K=2\), formed \(T\) (Hermitian, size 5) and \(G\) (real
symmetric, size 5), and checked whether
\(\lambda_{\min}(T)\ge -10^{-9}\) iff \(\lambda_{\min}(G)\ge
-10^{-9}\). **Zero mismatches**. The constraints are numerically
indistinguishable.

### Where the earlier verdict went wrong

* §3 ("Why \(G\) is a genuinely new constraint") contrasts \(G\)
  with the *moment lift* \(M\) (a free SDP block of v45), not with
  the *Hermitian Toeplitz* \(T\) (v46). That contrast is fine,
  but it is the wrong comparison: the question is whether \(G\)
  tightens v50 beyond what v50's own v46 block already enforces.
  It does not.
* The phrase "\(G\) involves \(a_{i+j}\) … whereas \(T\) only uses
  \(a_{|i-j|}\)" conflates **indexing** with the **set of moments
  used**. Both \(G\) and \(T^{(2K+1)}\) are built from the same
  moments \((a_k, b_k)_{k=0}^{2K}\); they just arrange them in
  different matrix patterns that turn out to be related by
  congruence.
* This is the same failure mode as Continuous Ramsey/v36:
  "unitary (here: non-unitary invertible) similarity of an
  existing PSD cut adds zero SDP information." (See
  `memory/feedback_verifier_rigor.md`.)

### What v50 would need to be a real tightening

Drop v46's size-\(2K+1\) Toeplitz from v50 and keep only v45's
size-\(K+1\) Toeplitz + the Gram \(G\) at size \(2K+1\). In that
setup, \(G\succeq 0\) would be a strict tightening of v45 (it's
equivalent to v46). But as currently written v50 = v45 + v46's
Toeplitz + \(G\), which is v45 + v46 + (tautology).

### Bottom line

The math of §3 of the earlier verdict (Gram entries, Gram PSD
proof) is correct, and the constraint itself is rigorously
derived. But the claim "\(G\succeq 0\) adds a new PSD cut on top
of v46" is wrong. If v50 ran numerically as described, it should
have given \(\Omega^{(50)}=\Omega^{(46)}\) to solver tolerance.
The project should confirm this by comparing numerics before
building further iterations on top of v50 as if it were a
non-trivial tightening.

(v51's localized Bochner with \(h_1=\cos(2\pi x)\) is *genuinely*
new relative to v46 — it uses the support restriction, which
neither v46 nor \(G\) does. That verdict stands.)

**Revised final verdict: INVALID AS A TIGHTENING (tautological
with v46). The underlying Gram-matrix math is correct; the claim
of novelty is wrong.**
