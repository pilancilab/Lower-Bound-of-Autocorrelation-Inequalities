# v48 — Higher-degree Fejér–Riesz via extended moment matrix

## Goal

In v44/v45 the Fejér–Riesz constraint acts on the degree-$K$ trig
polynomial $\Omega - \sigma_K^{f*f}(t)$.  The asymptote at $K=96$
appears to be around $\Omega \sim 1.294$ (extrapolating from the
v44 sweep: 1.289 @ N=1500, 1.290 @ N=2000, 1.291 @ N=3000, ...).
To break through this, we must increase the Fejér–Riesz **degree**.

The Cesàro error $\sup (f*f) - \sup \sigma_K^{f*f}$ is $O(1/K)$ for
Lipschitz $(f*f)$.  Doubling $K$ should halve the error.

## Construction

Define $K_{\text{big}} := 2K$.  Introduce:
- Extended Fourier variables $a_{K+1}, \dots, a_{K_{\text{big}}}$
  and $b_{K+1}, \dots, b_{K_{\text{big}}}$, with cell bounds as usual.
- Extended moment matrix $M_{\text{big}}$ of size
  $(2K_{\text{big}}+1) \times (2K_{\text{big}}+1) = (4K+1) \times (4K+1)$,
  lifting on $(1, a_1..K_{\text{big}}, b_1..K_{\text{big}})$.
- Extended Fejér–Riesz matrix $Q_{\text{big}}$ of size
  $(K_{\text{big}}+1) \times (K_{\text{big}}+1)$, Hermitian PSD,
  for the degree-$K_{\text{big}}$ trig polynomial.

The new Fejér–Riesz constraint:
$$
\Omega - \sigma_{K_{\text{big}}}^{f*f}(t) \ge 0 \quad \forall\,t,
$$
where $\sigma_{K_{\text{big}}}^{f*f}$ is the degree-$K_{\text{big}}$
Fejér mean with weights $w_k^{\text{big}} = 1 - k/(K_{\text{big}}+1)$,
encoded via $Q_{\text{big}}$ and $M_{\text{big}}[k,k]$,
$M_{\text{big}}[K_{\text{big}}+k, K_{\text{big}}+k]$,
$M_{\text{big}}[k, K_{\text{big}}+k]$ for $k = 1, \dots, K_{\text{big}}$.

The original $M$ and $Q$ (at degree $K$) can be **deleted** — they
are subsumed by $M_{\text{big}}$ and $Q_{\text{big}}$.  So v48 is
effectively v44 with $K \to 2K$.

## Why not just run v44 at $K = 2K$?

Running v44 directly at $K=192$ at $N=3000$ would have:
- Moment matrix size $385 \times 385$ (vs v44's $193 \times 193$)
- Fejér–Riesz $Q$ size $193 \times 193$ (vs v44's $97 \times 97$)

Solve cost scales as (matrix dim)$^3$, so roughly $8\times$ longer
per run.  At v44 (3000, 96) runtime of ~40 min, v48 (3000, 192) is
roughly **5 hours**.

**Compromise:** v48 = v44 at $K=128$ or $K=160$ at moderate $N$.

## Validity

Identical to v44.  The extended construction is a direct instance
of v44 with $K$ replaced by $K_{\text{big}}$.  All v44 validity
arguments carry over.

## What this is NOT

This is **not** a tightening of v44 in the sense of adding new
constraints; it is just v44 with larger $K$.  The reason for a
separate version is to emphasize the **parameter trajectory** and
to track whether the asymptote breaks the plateau.

## Expected results

At fixed moderate $N$ (e.g., $N = 500$):

| $K$ | est. solve time | expected $\Omega$ |
|-----|----------------:|------------------:|
| 96  | $\sim 15$ min   | $1.281$           |
| 128 | $\sim 40$ min   | $1.284$ ?         |
| 160 | $\sim 100$ min  | $1.286$ ?         |

## Status

Running.  Will update results table after the sweep.
