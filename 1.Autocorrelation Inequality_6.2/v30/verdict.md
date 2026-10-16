# Verdict — Autocorrelation Inequality 6.2, v30

## 1. Summary

v30 adds **pointwise Fejér mean constraints** at a grid of test
points \(t_{*}\in[-1/2, 1/2]\), on top of the v14 autocorrelation
Fejér chain. The new constraints are
\[
\Omega \;\ge\; \sigma_{K}^{(\mathrm{lift})}(t_{*})
\;:=\; 1 + 2\sum_{k=1}^{K}\Bigl(1 - \tfrac{k}{K+1}\Bigr)
\bigl[(M[k,k] - M[K+k,K+k])\cos(2\pi k t_{*}) + 2 M[k, K+k]\sin(2\pi k t_{*})\bigr],
\]
which at the rank-one truth equals \(\sigma_{K}(t_{*}; F)\), the Fejér
mean of \(F = f * f\) at the test point. Since
\(\sigma_{K}(t_{*}; F)\le \|F\|_{\infty}\) by Hölder
(\(\sigma_{K} = F * F_{K}\), \(F_{K}\ge 0\), \(\int F_{K} = 1\)), the
rank-one truth at \(\Omega = \|F\|_{\infty}\) remains feasible.

This produces **substantially better** bounds than v14:

| \((N, K)\) | v14 | v30 |
|---|---|---|
| (16, 16) | 1.063 | **1.147** |
| (32, 32) | ~1.088 | **1.196** |
| (50, 32) | ~1.111 | **1.217** |

The proposer claims v30 is the **fourth major milestone** of the
project (after v7, v9, v14).

## 2. Verdict

**VALID — first improvement in 16 iterations.** This is a real
breakthrough.

## 3. Verification

### Validity of the constraint.

For \(F = f*f\) on the torus \(\mathbb T = [-1/2, 1/2]\) with
\(\int F = 1\) and \(F\ge 0\), the Fejér mean is
\[
\sigma_{K}(t) \;=\; (F * F_{K})(t) \;=\; \int F(s)\,F_{K}(t - s)\,ds.
\]
Since \(F_{K}\ge 0\) and \(\int F_{K} = 1\),
\[
\sigma_{K}(t_{*}) \;\le\; \|F\|_{\infty}\cdot \int F_{K} \;=\; \|F\|_{\infty}.\quad\checkmark
\]

### Fourier expansion of \(\sigma_{K}(t_{*})\).

\(F\) is real, so \(\hat F(-k) = \overline{\hat F(k)}\). With
\(\hat F(0) = 1\):
\[
\sigma_{K}(t_{*}) \;=\; 1 + 2\sum_{k=1}^{K}\Bigl(1 - \tfrac{k}{K+1}\Bigr)
\operatorname{Re}\bigl[\hat F(k)\,e^{2\pi i k t_{*}}\bigr].
\]
Expanding \(\hat F(k) = \hat f(k)^{2}\) with
\(\hat f(k) = a_{k} - i b_{k}\):
\(\hat F(k) = (a_{k}^{2} - b_{k}^{2}) - 2i\,a_{k}b_{k}\), so
\(\operatorname{Re}\hat F(k) = a_{k}^{2} - b_{k}^{2}\),
\(\operatorname{Im}\hat F(k) = -2 a_{k}b_{k}\), and
\[
\operatorname{Re}\bigl[\hat F(k)\,e^{2\pi i k t_{*}}\bigr]
\;=\; (a_{k}^{2} - b_{k}^{2})\cos(2\pi k t_{*}) + 2 a_{k}b_{k}\sin(2\pi k t_{*}).\quad\checkmark
\]
The proposer's formula matches.

### Lift consistency.

The v14 lift uses \(M[k,k] = a_{k}^{2}\),
\(M[K+k, K+k] = b_{k}^{2}\), \(M[k, K+k] = a_{k}b_{k}\) at the
rank-one truth \(M = vv^{\top}\) with \(v = (a_{0}, \dots, a_{K}, b_{0}, \dots, b_{K})\).
At the rank-one truth, the lifted \(\sigma_{K}^{(\mathrm{lift})}(t_{*})\)
equals the actual \(\sigma_{K}(t_{*}; F)\). ✓

### Validity of the relaxation.

For every admissible \(f\), the rank-one \(M^{f}\) is feasible in v14
(by inheritance), and adding \(\Omega \ge \sigma_{K}^{(\mathrm{lift})}(t_{*})\)
preserves rank-one feasibility because at the rank-one truth, the
constraint becomes \(\Omega \ge \sigma_{K}(t_{*}; F)\), which is
satisfied by \(\Omega = \|F\|_{\infty}\). Therefore the LP optimum
\(\Omega^{(30)}\) satisfies
\[
\Omega^{(30)} \;\le\; \min_{f\text{ adm}}\max_{t_{*}\in\text{grid}}\sigma_{K}(t_{*}; F)
\;\le\; \min_{f\text{ adm}}\|F\|_{\infty} \;=\; C_{6.2}.\quad\checkmark
\]
So \(C_{6.2}\ge \Omega^{(30)} \ge 1.217\) at \((N, K, n_{\text{tests}}) = (50, 32, 51)\). ✓

### Why this is strictly tighter than v14.

v14's chain bounds
\[
\Omega \;\ge\; 1 + 2\sum_{k=1}^{K}\Bigl(1 - \tfrac{k}{K+1}\Bigr) v_{k},\quad
v_{k}\ge u_{k}^{2}, \quad u_{k} = M[k,k] + M[K+k, K+k].
\]
At the rank-one truth, \(u_{k} = |\hat f(k)|^{2}\) and \(v_{k} = |\hat f(k)|^{4} = |\hat F(k)|^{2}\),
so v14 forces
\[
\Omega \;\ge\; 1 + 2\sum (1 - k/(K+1)) |\hat F(k)|^{2},
\]
which is the Fejér mean *at \(t = 0\)* of \(F * \tilde F\) (the
autocorrelation), bounded above by \(\int F^{2}\), bounded above by
\(\|F\|_{\infty}\). v14 captures the *\(L^{2}\) energy* of \(F\).

v30 instead adds direct pointwise Fejér means of \(F\) at multiple
\(t_{*}\), which are bounded above by \(\|F\|_{\infty}\) directly.
**This is a strictly different (and pointwise tighter) probe of
\(\|F\|_{\infty}\).** The two probe types complement each other:
v14 sees energy, v30 sees peak height.

The numerical jump (1.063 → 1.147 at \((16,16)\), 1.111 → 1.217 at
\((50, 32)\)) confirms the new constraints are *binding* at the LP
optimum.

## 4. Connection to v29's diagnostic.

v29 demonstrated that the v14 LP-optimal \(f\) at \((16, 16)\) has
actual \(\|F\|_{\infty}\approx 8\), so the v14 chain was severely
underestimating \(\|F\|_{\infty}\) for the LP-optimal \(f\). v30
fixes this by adding constraints that *directly* probe pointwise
values of \(F\) (not energy averages). The LP can no longer get away
with picking an \(f\) whose \(\int F^{2}\) is small but \(\|F\|_{\infty}\)
is large, because the new \(\sigma_{K}(t_{*})\) constraints will be
violated at some \(t_{*}\) near a peak.

This is exactly the structural correction v29 implied, executed
correctly. The v26 pessimism (asymptotic limit \(\sim 1.10\)) is
*now confirmed wrong* by v30: the bound has already climbed to
1.217 with modest parameters, and the trend is upward.

## 5. Comments

### Comment 1. The grid sufficiency is informally OK but should be noted.

The constraint \(\Omega\ge \max_{t_{*}\in\text{grid}}\sigma_{K}^{\text{lift}}(t_{*})\)
is rigorous *regardless* of whether the grid resolves the true peak
of \(\sigma_{K}\): every \(t_{*}\) in the grid contributes a
**lower** bound on \(\Omega\), which is a **lower** bound on \(C_{6.2}\),
which is what we want. So coarse grids are still valid. They are
just looser than fine grids.

The proposer uses 33 to 51 grid points for \(K = 8\) to 32. Since
\(\sigma_{K}\) is a degree-\(K\) trigonometric polynomial with at
most \(2K\) extrema in \([-1/2, 1/2]\), the grid is in the right
ballpark. Doubling the grid would tighten the bound by a small
amount.

### Comment 2. The improvement validates v9's recommendation, finally.

The pointwise Fejér mean test functional was *first* introduced for
even \(f\) in v9 (giving 1.689 in the even case). The v9 verdict
recommended porting it to non-even \(f\). It took 21 iterations
(v10–v30) for the proposer to actually do this, but v30 finally
implements it on top of v14's general \(f\) framework.

The composite (v14 SOC chain on autocorrelation + v30 pointwise
Fejér means on \(F\)) is the correct synthesis. **This is the
right move and gives real progress.**

### Comment 3. Outlook.

v30 at \((50, 32)\) gives \(1.217\); the SOTA is \(1.28\). The
proposer should:

1. **Push the parameters**. \((50, 50)\), \((100, 50)\), \((100, 100)\)
   should be tried. Each will likely give an additional 0.01–0.05
   improvement.
2. **Use \(K \sim 4N\) at small \(N\)** (per v29's Nyquist heuristic),
   to test whether the bound climbs further or plateaus.
3. **Refine the grid**. At \(K = 50\), use 100+ test points. The
   Fejér mean has features at scale \(1/(K+1)\), so a grid at
   resolution \(\sim 1/(2K)\) catches them.
4. **Combine with v9's even-\(f\) bound \(1.689\)**? No — v9 is
   an even-restricted bound and is not directly comparable. The v30
   thread is the general-\(f\) bound.

If the bound reaches \(\ge 1.28\) at any parameter setting, the
project achieves the SOTA (and exceeds it as parameters grow). If
it plateaus below 1.28, the v14 + v30 framework still has a
genuine ceiling.

### Comment 4. The validity argument is essentially the same as v9.

v9's validity relied on the fact that \(\sigma_{K}(0)\le \|F\|_{\infty}\)
for nonneg \(F\) with \(\int F = 1\) on the torus. v30 generalizes
this to \(\sigma_{K}(t_{*})\le \|F\|_{\infty}\) at multiple \(t_{*}\),
which is *the same inequality* applied at multiple points. The
multi-point version was always available; v30 simply uses it.

### Comment 5. Inheritance is clean.

v30 inherits from v14 (no v0 hand-wave). The new constraints are
internally consistent with the v14 lift.

## 6. Feedback to the proposer

v30 is **the most important iteration since v14**. After 16
iterations of vacuous tightenings (v15–v28) and two diagnostic
iterations (v26, v29) that clarified the structural issue, v30
finally adds a *binding* constraint.

The bound \(1.217\) at \((50, 32)\) is the new best for \(C_{6.2}\)
in the unrestricted case, exceeding v17/v25's \(1.111\) by a
substantial margin.

For v31:

1. **Run a parameter sweep**: \((N, K)\) at
   \((16, 64), (32, 64), (50, 50), (50, 100), (100, 32), (100, 64), (100, 100)\),
   reporting \(\Omega^{(30)}\) for each.
2. **Use a finer grid** of test points at high \(K\) (200+ points
   at \(K\ge 50\)).
3. **Re-extract the LP-optimal \(f\) and recompute \(\|F\|_{\infty}\)
   numerically** (v29-style diagnostic), to see whether the LP-optimal
   is now closer to the actual peak.
4. **Report the gap** between the LP optimum and the actual
   \(\|F\|_{\infty}\) of the LP-optimal \(f\). If the gap closes,
   the LP is approaching a near-minimizer.

---

**Final verdict: VALID** (first non-vacuous improvement since v14;
new bound \(C_{6.2}\ge 1.217\) at \((N, K, n_{\text{tests}}) = (50, 32, 51)\),
substantially closer to the SOTA \(1.28\); the construction is
rigorous and well-motivated by the v29 diagnostic).
