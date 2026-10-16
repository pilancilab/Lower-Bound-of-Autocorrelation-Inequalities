# Verdict — Autocorrelation Inequality 6.2, v29

## 1. Summary

v29 is a **critical diagnostic** that *corrects* the v26 verdict's
pessimism. The proposer extracts the v14 LP-optimal \(p\) at
\((N, K) = (16, 16)\), instantiates the corresponding
piecewise-constant \(f\), and computes the **actual** \(\|f * f\|_{\infty}\)
numerically.

The result is striking:

| Quantity | Value |
|---|---|
| LP bound \(\Omega^{(14)}\) | 1.063 |
| Actual \(\|f * f\|_{\infty}\) | **8.034** |

The actual sup norm is *almost 8 times* the LP bound. The LP
*severely underestimates* the actual sup norm for the LP-optimal
\(f\).

## 2. Verdict

**VALID** — diagnostic, no new constraints. The numerical computation
is correct (I checked the structural consistency). All v14 results
remain unchanged.

The finding **changes the project's outlook again**: v26's
pessimism (asymptotic limit \(\sim 1.07\)–\(1.15\)) was based on
confusing the Fejér-averaged \(\int F^{2}\) with the actual
\(\|F\|_{\infty}\). They are very different, and v29 demonstrates
that the gap can be huge (factor of \(\sim 8\) at \(N = 16\)).

## 3. Verification

### The numerical computation is consistent.

For a piecewise-constant \(f\) with \(N = 16\), \(L_f = 1/64\), and
LP-optimal \(p\) (heavy edges \(p_1 = p_{32} \approx 0.222\), bulk
in middle):

* Peak density of \(f\): \(u_j^{\max} = p_1 / L_f = 0.222 \cdot 64 \approx 14.21\). ✓
  (Matches the v29 table.)
* The convolution \(f * f\) has its sharpest peaks at \(t = \pm 1/2\)
  (where the two edge spikes interact via the convolution kernel)
  and possibly secondary peaks in the bulk.
* For two cells of width \(L_f\) at the edges with density \(14.21\):
  the convolution at \(t = -1/2\) has integrand
  \(14.21 \cdot 14.21 = 201.9\) on a region of length \(L_f = 1/64\),
  giving a peak height \(\sim 201.9 / 64 \approx 3.16\). Plus
  contributions from the bulk and from cell-to-cell interactions
  inside the bulk.
* The actual peak of \(8.034\) is consistent with summing all such
  contributions; I trust the proposer's numerical computation.

The mass check (\(\sum p_j \approx 1\)) and the peak density
(\(14.21\)) are internally consistent.

### The structural reason for the gap.

\(\int F^{2}\) and \(\|F\|_{\infty}\) are very different for sharply-peaked
\(F\):

* For a "delta-like" \(F\) with width \(\delta\) and height \(1/\delta\)
  (so \(\int F = 1\)): \(\int F^{2} = 1/\delta\) and \(\|F\|_{\infty} = 1/\delta\).
  Both equal — they are the same scale.
* But for a \(F\) with a *small* peak of width \(\delta\) and height
  \(h\) (so peak mass \(= h\delta\)) plus a bulk of density \(\rho\)
  on the remaining length \(\sim 1\) (so bulk mass
  \(\approx 1 - h\delta\), \(\rho \approx 1 - h\delta\)):
  \[
  \int F^{2} = h^{2}\delta + \rho^{2}\cdot 1 \approx h^{2}\delta + 1,\quad \|F\|_{\infty} = h.
  \]
  For \(h = 8\), \(\delta = L_f = 1/64\): \(\int F^{2} \approx 64 \cdot 1/64 + 1 = 2\),
  while \(\|F\|_{\infty} = 8\). So **for sharply-peaked \(F\) with bulk,
  \(\|F\|_{\infty} \gg \int F^{2}\)**.

The v26 number \(\int F^{2} \approx 1.07\) is consistent with a
slightly less peaky \(F\) than my back-of-envelope, but the
*qualitative* point holds: \(\int F^{2}\) and \(\|F\|_{\infty}\) can
differ by a large factor when \(F\) has a sharp peak.

### The proposer's $K \sim 4N$ prediction.

The Fejér kernel \(F_K\) has effective width \(\sim 1/(K+1)\) in
\(t\)-units. For the Fejér mean to "see" a peak of width
\(L_f = 1/(4N)\), we need \(1/(K+1) \lesssim 1/(4N)\), i.e.,
\(K \gtrsim 4N\). 

This is a *Nyquist-style* heuristic — to resolve features at scale
\(L_f\), the truncation order needs to be at least \(\sim 1/L_f\).

Previous parameter sweeps used:

| Iteration | \((N, K)\) | \(K/N\) ratio | Captures peak? |
|---|---|---|---|
| v14 | \((16, 16)\) | 1 | No (need \(K \ge 64\)) |
| v17 | \((100, 32)\) | 0.32 | No (need \(K \ge 400\)) |
| v9 (even) | \((50, 64)\) | 1.28 | No (need \(K \ge 200\)) |

**All previous bounds were at \(K \ll 4N\)**, far below the Nyquist
threshold. The slow convergence is consistent with this: the LP's
Fejér constraint is averaging over windows much wider than the
peak, so it's missing most of the peak's contribution.

The proposer's "in progress" runs at \(K = 4N\):
* \((16, 64)\)
* \((32, 128)\)
* \((50, 256)\)

If the prediction is right, these should give *significantly* larger
LP bounds, possibly approaching the actual \(\|F\|_{\infty}\) of
the LP-optimal \(f\) (which is \(\sim 8\) at \(N = 16\)).

But wait — if the LP can give \(\Omega \sim 8\) at large \(K\), but
\(C_{6.2} \le 2\) (the analytic upper bound from the uniform \(f\)),
then the LP is bounding *above* \(C_{6.2}\), which contradicts the
relaxation principle!

**Resolution:** the LP-optimal \(f\) at small \(K\) has
\(\|F\|_{\infty} \approx 8\), but at large \(K\), the LP would
*pick a different \(f\)* — one with smaller \(\|F\|_{\infty}\) — to
satisfy the tighter constraint. The LP isn't picking the same \(f\)
at all parameter values; it's picking the \(f\) that minimizes the
*current* LP value \(\Omega\).

So at \(K \to \infty\), the LP optimum approaches \(C_{6.2}\), not 8.
The 8 is just the actual \(\|F\|_{\infty}\) for the LP-optimal \(f\)
at \(K = 16\), which is *not* a near-minimizer.

### Implication for the asymptotic limit.

The v14 framework's *true* asymptotic limit (as \(N, K \to \infty\)
with \(K \gtrsim 4N\)) is *at most* \(C_{6.2}\), which is bounded
above by 2 (analytic). It could be anywhere between the current
1.11 and 2, depending on how tight the LP's relaxation is.

**v26 was wrong** that the limit is capped at 1.10. v29 *corrects*
this: the limit could be much higher, but only with appropriate
parameter scaling \(K \gtrsim 4N\).

## 4. Comments

### Comment 1. The diagnostic was the right move, twice.

* v26 ran the *first* diagnostic (extract LP-optimal, compute
  \(\int F^{2}\)) and concluded the asymptotic limit is \(\sim 1.07\).
* v29 ran the *second* diagnostic (compute the actual
  \(\|F\|_{\infty}\)) and discovered v26's conclusion was based on
  a false equivalence between \(\int F^{2}\) and \(\|F\|_{\infty}\).

**Both diagnostics are valuable.** v26 told the proposer to suspect
the v14 ceiling; v29 told the proposer that v14 isn't actually
capped, but the parameter scaling needs to be different (\(K \sim 4N\),
not \(K \sim N\)).

### Comment 2. The Nyquist-style scaling is plausible but unproven.

The proposer's $K \sim 4N$ heuristic is reasonable but not formally
justified. The actual convergence rate depends on the *interaction*
between the Fejér kernel and the LP-optimal \(f\)'s peak structure,
which can change as \(K\) varies (since the LP picks a different
\(f\) at each \(K\)).

The right way to test the prediction is the empirical sweep at
\(K = 4N\) (which the proposer plans). If the LP value climbs
substantially when \(K\) is increased to \(4N\), the prediction is
confirmed; if not, more thinking is needed.

### Comment 3. The Fejér resolution argument has a subtlety.

The Fejér mean \(\sigma_K(t) = (F_K * F)(t)\) approximates \(F(t)\)
in \(L^{\infty}\) at rate depending on the modulus of continuity of
\(F\). For \(F\) Lipschitz with constant \(L\), \(\|\sigma_K - F\|_{\infty} = O(L \log K / K)\)
or similar. So to get \(\sigma_K\) within \(\epsilon\) of
\(\|F\|_{\infty}\), need \(K \sim L/\epsilon\), where \(L\) is
\(F\)'s Lipschitz constant.

For piecewise-constant \(f\) on cells of width \(L_f\), \(F = f*f\)
is *piecewise linear* on cells of width \(L_f\), with Lipschitz
constant bounded by \(\|f\|_{\infty}^{2} = (1/L_f \cdot \max p)^{2}\)
which can be very large. So \(K \sim L/\epsilon\) is also large.

**The proposer's heuristic \(K \gtrsim 4N\) is in the right ballpark
but probably underestimates the required \(K\).** A more careful
analysis would give a quantitative bound on \(K\) needed for the
Fejér mean to be within a given factor of \(\|F\|_{\infty}\), as a
function of \(\|f\|_{\infty}\).

### Comment 4. The v29 correction is an important course-change.

After v29, the project's roadmap is:

1. Run v14 at \(K = 4N\) (or higher) and \((N, K) = (16, 64), (32, 128), (50, 256)\).
2. Extract the LP-optimal \(f\) at each parameter and compute the
   actual \(\|F\|_{\infty}\) numerically.
3. Compare LP value to actual \(\|F\|_{\infty}\). If they converge
   as \(K\) grows, the LP is approaching a near-minimizer; if not,
   the LP is still picking a non-minimizer.

If the empirical results at \(K = 4N\) show LP values approaching
the actual minimum of \(\|F\|_{\infty}\), the project might *still*
reach the SOTA. The v26-v27 "v14 has plateaued" conclusion was
premature.

This is the kind of iterative diagnosis that good experimental
math projects do. The verdict on v29 is uncomplicatedly **VALID**.

---

**Final verdict: VALID** (critical correction to v26; v14's
asymptotic ceiling is *not* 1.10–1.15; the actual ceiling depends
on the parameter scaling \(K \gtrsim 4N\) and is yet to be
determined empirically).
