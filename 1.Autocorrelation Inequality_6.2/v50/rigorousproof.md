# v50 — Gram-matrix constraint in the real trig basis

## Goal

Add a genuinely new PSD constraint to the v45 construction:
$$
G \;:=\; \Bigl[\,\int f(x)\,g_i(x)\,g_j(x)\,dx\,\Bigr]_{i,j=0}^{2K} \;\succeq\; 0,
$$
where $(g_i)_{i=0}^{2K}$ is the **real** trigonometric basis
$$
g_0 = 1,\qquad
g_k = \cos(2\pi k x)\ (1 \le k \le K),\qquad
g_{K+k} = \sin(2\pi k x)\ (1 \le k \le K).
$$
$G$ is the **Gram matrix** of this basis in the $f$-weighted inner
product $\langle \phi, \psi\rangle_f := \int f\,\phi\,\psi$, which
is positive semi-definite whenever $f \ge 0$ — automatic regardless
of the SDP relaxation slack.

Crucially, **the entries of $G$ are LINEAR in the Fourier moments**
$(a_k, b_k)$ for $k = 0, 1, \dots, 2K$, via the standard trig
product-to-sum identities.  So $G$ becomes a new PSD constraint on
the extended moment sequence, without any additional lifting.

## Entries of $G$

Using the identities
$$
\cos\alpha \cos\beta = \tfrac{1}{2}(\cos(\alpha-\beta) + \cos(\alpha+\beta)),
$$
$$
\sin\alpha \sin\beta = \tfrac{1}{2}(\cos(\alpha-\beta) - \cos(\alpha+\beta)),
$$
$$
\cos\alpha \sin\beta = \tfrac{1}{2}(\sin(\alpha+\beta) - \sin(\alpha-\beta)),
$$
we get:

- $G[0,0] = \int f = 1$
- $G[0, i] = G[i, 0] = \int f \cos(2\pi i x) = a_i$ for $1 \le i \le K$
- $G[0, K+j] = G[K+j, 0] = \int f \sin(2\pi j x) = b_j$ for $1 \le j \le K$
- For $1 \le i, j \le K$ (cos-cos block):
  $$
  G[i, j] \;=\; \int f \cos(2\pi i x)\cos(2\pi j x)\,dx
  \;=\; \tfrac{1}{2}\bigl(a_{|i-j|} + a_{i+j}\bigr).
  $$
- For $1 \le i, j \le K$ (sin-sin block):
  $$
  G[K+i, K+j] \;=\; \int f \sin(2\pi i x)\sin(2\pi j x)\,dx
  \;=\; \tfrac{1}{2}\bigl(a_{|i-j|} - a_{i+j}\bigr).
  $$
- For $1 \le i, j \le K$ (cos-sin block):
  $$
  G[i, K+j] \;=\; \int f \cos(2\pi i x)\sin(2\pi j x)\,dx
  \;=\; \tfrac{1}{2}\bigl(b_{i+j} + \mathrm{sign}(j-i)\,b_{|j-i|}\bigr),
  $$
  where $\mathrm{sign}(0) := 0$ (so the diagonal $i = j$ gives $G[i, K+i] = b_{2i}/2$).

All entries use $a_k, b_k$ for $k \in \{0, 1, \dots, 2K\}$.

## Validity

**Claim.** $G \succeq 0$ holds for every rank-1 feasible point,
i.e., for every real $f \ge 0$ with finite mass on
$[-1/4, 1/4]$.

**Proof.** For any real vector $v = (v_0, v_1, \dots, v_{2K})$,
$$
v^\top G v \;=\; \sum_{i, j} v_i\,v_j\,\int f\,g_i\,g_j\,dx
\;=\; \int f(x)\, \Bigl(\sum_i v_i\,g_i(x)\Bigr)^2\,dx \;\ge\; 0,
$$
since $f \ge 0$ and the integrand is a nonneg function.  Hence
$G \succeq 0$.  $\square$

**Corollary.** $G \succeq 0$ is a valid relaxation constraint to add
to the SDP.  Every rank-1 feasible point satisfies it, so adding it
excludes no genuine feasible $(a, b, M, Q)$.

## Is this constraint new?

Compared to the existing v45 constraints:

1. **$M \succeq 0$ (moment lift)** acts on the variables
   $(1, a_1, \dots, a_K, b_1, \dots, b_K)$ via the PSD relaxation of
   their rank-1 outer product.  Entries $M[i,j]$ are SDP variables.
2. **Bochner–Toeplitz $T \succeq 0$** acts on the Hermitian Toeplitz
   with first row $c_k = a_{|k|} - i\,\mathrm{sgn}(k)\,b_{|k|}$.
   Entries of $T$ are linear in $(a_k, b_k)$ up to $k = K$ (or $2K$
   in v46).
3. **$G \succeq 0$** (new): acts on the **real** trig basis with
   entries linear in $(a_k, b_k)$ for $k = 0, \dots, 2K$, using the
   product-to-sum identities.

The $G$ matrix is **different** from $T$ because the real trig
basis $\{1, \cos, \sin\}$ produces entries involving both $a_{|i-j|}$
and $a_{i+j}$ (from the product-to-sum), whereas the Hermitian
Toeplitz $T$ only uses $a_{|i-j|}$ (from exponentials
$e^{\pm 2\pi i (i-j) x}$).

Thus $G \succeq 0$ encodes positivity information about **higher-
frequency** Fourier moments $a_{i+j}$ for $2 \le i+j \le 2K$, which
is not captured by the lower-frequency Toeplitz structure alone.

At the **rank-1 limit** (true $f$), $G$ coincides with the
moment-matrix lift
$M_{\text{rank-1}}[i,j] = g_i(x_0)\,g_j(x_0)$ for a single-atom $f$;
the identity $g_i g_j = \frac{1}{2}(g_{|i-j|} + g_{i+j})$ is exactly
the rank-1 reduction via the trig product rule.

## Relationship to the C–S constraints of v47

v47 added $M[k,k] \le G[k,k] = (1 + a_{2k})/2$ and its $b$-analog.
These are **diagonal** entries of $G$ and are subsumed by the full
PSD constraint $G \succeq 0$.

The v47 diagonal bounds were **not binding** (the SDP pushes $M[k,k]$
down, not up).  But the **full** $G \succeq 0$ constraint may be
different, because it constrains the off-diagonal structure of $G$
— a constraint on combinations of $(a_k, b_k, a_{i+j}, b_{i+j})$
that the SDP cannot avoid.

## Construction

v50 = v45 + (aux variables $a_{K+1..2K}$, $b_{K+1..2K}$ with cell
bounds, as in v46) + ($G \succeq 0$).

The extended Bochner–Toeplitz of v46 can optionally be kept or
dropped.  v50 keeps it, since it's cheap.

## Status

**Awaiting theorist verdict.**  The new constraint is a standard
Gram-matrix PSD of a basis in an $f$-weighted inner product —
textbook-level validity.  The novelty is its specific form in the
real trig basis with the product-to-sum coupling between $k$ and
$2K$.

## Expected gain

Hard to predict.  If the SDP solution's $M$ entries are already
forced into a rank-1-like shape by the other constraints, adding
$G \succeq 0$ may not help.  But if the relaxation is loose in the
cos-sin block (where $G$ has entries $b_{2i}/2$ that $T$ does not
directly constrain), there should be some tightening.
