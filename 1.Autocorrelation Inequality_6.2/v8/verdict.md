# Verdict — Autocorrelation Inequality 6.2, v8

## 1. Summary of the proof

v8 attempts the F-side / Fourier-side relaxation suggested in the v7
verdict. It introduces:

* f-side cell averages \(p_j\) (inherited from v7),
* F-side cell averages \(W_m\),
* Fourier coefficients \(a_k\) of \(f\) (even \(f\) only),
* a PSD lift \(M\) of the moment vector \((1, a_1, \ldots, a_K)\),
* a Fourier link \(\hat F(k) = \hat f(k)^{2}\) giving \(M[k,k] \approx A_k(F)\),
* and — critically — a **naive Parseval lower bound**
  \(\sum_{k=1}^{K} M[k,k] \ge 1/2\) intended to capture
  \(\sum_{k\ge 1} a_k^{2} = (\int f^{2} - 1)/2 \ge 1/2\).

The proposer then **catches their own mistake** before submitting:
they verify that the uniform admissible \(f = 2\cdot\mathbf 1_{[-1/4, 1/4]}\)
violates the truncated Parseval constraint (the truncated sum
\(\sum_{k=1}^{4} a_k^{2} \approx 0.45 < 0.5\)), and explicitly mark v8
as INVALID with a v9 plan to use the **Fejér kernel test functional**
instead.

## 2. Verdict

**INVALID** (self-INVALID, correctly diagnosed by the proposer).

The proposer is right. The naive Parseval lower bound at finite \(K\)
is *not* a valid LP constraint for admissible \(f\), and the proof
correctly identifies the admissible counterexample (uniform \(f\)),
explains the structural reason (no a-priori bound on \(\int f^{2}\)
in the admissible set), and proposes the correct fix (Fejér kernel).

This is *exactly* the kind of self-checking the verifier instructions
ask for. **v8 is honestly self-INVALID, and that is the correct
verdict.**

## 3. Verification of the self-INVALID

I checked the proposer's calculation independently and confirm.

### The Fourier coefficients of uniform \(f\). ✓

For \(f \equiv 2\) on \([-1/4, 1/4]\) (so \(\int f = 1\)):
\[
a_k \;=\; \int_{-1/4}^{1/4} 2\cos(2\pi k x)\,dx
\;=\; 2\Big[\frac{\sin(2\pi k x)}{2\pi k}\Big]_{-1/4}^{1/4}
\;=\; \frac{2\sin(\pi k/2)}{\pi k}.
\]

For \(k\) even: \(\sin(\pi k/2) = 0\), so \(a_k = 0\).
For \(k\) odd: \(|a_k| = 2/(\pi k)\), so \(a_k^{2} = 4/(\pi^{2} k^{2})\).

### The truncated sum at \(K = 4\). ✓
\[
\sum_{k=1}^{4} a_k^{2}
\;=\; a_1^{2} + 0 + a_3^{2} + 0
\;=\; \frac{4}{\pi^{2}} + \frac{4}{9\pi^{2}}
\;=\; \frac{36 + 4}{9\pi^{2}}
\;=\; \frac{40}{9\pi^{2}}
\;\approx\; 0.4502.
\]
This is **strictly less than \(1/2 = 0.5\)**, so the LP constraint
\(\sum_{k=1}^{4} M[k,k] \ge 1/2\) **excludes** the rank-one lift from
the uniform admissible \(f\). Therefore the constraint is *not*
satisfied by every admissible \(f\), and v8's claim
\(\Omega^{(8)}\le C^{\text{even}}_{6.2}\) is not certifiable. ✓

### Sanity check: the *infinite* Parseval identity holds. ✓
\[
\sum_{k=1}^{\infty} a_k^{2} = \sum_{k\text{ odd}}\frac{4}{\pi^{2} k^{2}}
\;=\; \frac{4}{\pi^{2}}\sum_{k\text{ odd}}\frac{1}{k^{2}}
\;=\; \frac{4}{\pi^{2}}\cdot\frac{\pi^{2}}{8}
\;=\; \frac{1}{2},
\]
using the standard \(\sum_{k\text{ odd}} 1/k^{2} = \pi^{2}/8\). And
\(\int f^{2} = 4 \cdot 1/2 = 2\), giving
\((\int f^{2} - 1)/2 = 1/2\). ✓

So the infinite-sum identity is correct, but the truncation
\(\sum_{k=1}^{K}\) is *strictly less* than \(1/2\) for any finite
\(K\) (the tail \(\sum_{k > K}\) is positive but truncated to zero).
The naive truncation drops a positive term and so cannot give a
*lower* bound on the partial sum.

This is the *same* tail-bound problem I flagged in **6.2 v2 Issue 2**
and **6.2 v6 Comment 1**: any constraint of the form "the truncated
Fourier-side moment is bounded below" requires an explicit tail
control, and there is no a-priori upper bound on \(\int f^{2}\) in
the admissible set (so no a-priori bound on the tail).

### Why no a-priori bound on \(\int f^{2}\). ✓

The proposer correctly notes:
> \(\int f^{2}\) is unbounded above on the admissible set (e.g.,
> \(f = c\cdot\mathbf 1_{[-1/(4c), 1/(4c)]}\) for \(c\to\infty\) has
> \(\int f^{2}\to\infty\)).

I checked: \(\int f = c\cdot 2\cdot\frac{1}{4c} = 1/2\). Hmm, that's
not 1. Let me re-check: \(f = c\cdot\mathbf 1_{[-1/(4c), 1/(4c)]}\)
has \(\int f = c\cdot \frac{1}{2c} = 1/2\). To make this admissible
with \(\int f = 1\), use \(f = 2c\cdot\mathbf 1_{[-1/(4c), 1/(4c)]}\),
which has \(\int f = 2c\cdot\frac{1}{2c} = 1\) ✓ and
\(\int f^{2} = (2c)^{2}\cdot\frac{1}{2c} = 2c \to \infty\) as
\(c\to\infty\). The proposer's example is morally correct (just off
by a factor of 2 in the height); the conclusion that
\(\int f^{2}\) is unbounded on the admissible set is correct.

(Alternatively: a delta-like sequence
\(f_\epsilon = (1/\epsilon)\mathbf 1_{[-\epsilon/2, \epsilon/2]}\) on
\([-1/4, 1/4]\) has \(\int f_\epsilon = 1\) and
\(\int f_\epsilon^{2} = 1/\epsilon\to\infty\). Same conclusion.)

So the tail \(\sum_{k > K} a_k^{2}\) can be made arbitrarily large
by concentrating \(f\), and no LP-side variable can bound it.

## 4. Comments

### Comment 1. The self-check is exactly what the verifier instructions ask for.

The verifier instructions say:

> MATHEMATICAL RIGOR IS THE MOST IMPORTANT!!! Do not simply assume
> that you are doing an approximation that is good enough or do not
> propose things that may or may not be a mathematically valid lower
> or upper bound.

The earlier iterations (v0, v1, v2, v3, v4, v5, v6) violated this by
inheriting hand-waved chains without testing them on admissible
counterexamples. v8 *does* the test (against the uniform \(f\)),
*finds* the counterexample, and *honestly* marks itself INVALID. This
is the right scientific practice and is to be commended.

### Comment 2. The v9 Fejér kernel plan is the correct fix.

The proposer's v9 plan:

> v9 will:
> 1. Replace the naive Parseval with a Fejér kernel test functional
>    that gives a rigorous linear inequality without any tail issue.
> 2. Use the cleaner \(\int F_K(t) F(t)\,dt \le \Omega \cdot \int F_K(t)\,dt = \Omega\)
>    from the verifier's recommendation, where \(F_K\) is the (positive)
>    Fejér kernel of order \(K\).

is *exactly* what I recommended in the **6.2 v2 verdict Section 4**,
**6.2 v6 verdict Section 4**, and (in spirit) the **6.2 v7 verdict
Section 5**. The Fejér kernel
\(F_K(t) = \frac{1}{K+1}\Big(\frac{\sin((K+1)\pi t/T)}{\sin(\pi t/T)}\Big)^{2}\)
(or its discrete analog) is **nonnegative**, so for any admissible
\(F = f * f\),
\[
\int F_K(t)\,F(t)\,dt \;\le\; \Omega \cdot \int F_K(t)\,dt,
\]
with no tail correction (the Fejér kernel has finite Fourier support,
so the LHS is exactly a partial sum of the Fourier coefficients of
\(F\) — and there is no Gibbs overshoot to compensate for, because
\(F_K\) is nonneg and the inequality direction is correct).

This gives a **rigorous linear inequality** in the Fourier
coefficients of \(F\) (which, via the squared identity
\(\hat F(k) = \hat f(k)^{2}\), are the diagonal entries
\(M[k,k] = a_k^{2}\) of the moment matrix). And it requires *no
tail bound* on \(\int f^{2}\), because the Fejér mean already
truncates the Fourier expansion at \(K\) by construction.

I strongly recommend this for v9, and I am confident it will produce
a non-trivial constraint that *does* slot into the LP/SDP cleanly.

### Comment 3. The proof's "Numerical observation" is a small bonus.

Even ignoring the (P4) rigor failure, the proof reports that MOSEK
returns \(\Omega^{(8)} = 1\) — i.e., the constraint as currently
written *also* fails to bind numerically. This is a useful sanity
check: even if the constraint had been rigorous, it wouldn't have
moved the bound. The Fejér kernel approach should fix both issues
(rigor and binding) simultaneously.

### Comment 4. The inherited v0 band hand-wave issue is now moot.

v7 dropped the v0 band constraint and replaced it with the rigorous
(W) family. v8 inherits v7, *not* v0, so the inherited band-constraint
hand-wave is gone. v8's rigor failure is a *new* one (the (P4)
Parseval truncation), independent of the old issue.

## 5. Feedback to the proposer

**This is the right kind of failure.** Catching the bug *before*
submitting, identifying the structural reason (no a-priori bound on
\(\int f^{2}\)), and proposing the correct fix (Fejér kernel) is
exactly the workflow I'd want to see. The verdict on v8 is
**INVALID**, but the *project* is in much better shape than it would
be if v8 had been submitted without the self-check.

Two concrete suggestions for v9:

1. **Use the Fejér kernel test functional cleanly.** Concretely, for
   each test point \(t_*\in[-1/2, 1/2]\), the Fejér mean of \(F\) at
   \(t_*\) is
   \[
   \sigma_K(t_*) \;=\; \sum_{k=-K}^{K}\Big(1 - \frac{|k|}{K+1}\Big)\hat F(k)\,e^{2\pi i k t_*}.
   \]
   By the standard Fejér-mean theorem (which uses the *nonnegativity*
   of the Fejér kernel), \(\sigma_K(t_*) \ge 0\) and
   \(\sigma_K(t_*) \le \|F\|_{\infty} = \Omega\). And by the squared
   identity \(\hat F(k) = \hat f(k)^{2}\), \(\sigma_K(t_*)\) is a
   *bilinear* function of \(\hat f\), which lifts to a *linear*
   function of the moment matrix \(M\). So
   \(\sigma_K(t_*) \le \Omega\) is a linear constraint on \(M\).
   This is **rigorous** with **no tail correction**.

2. **Test the v9 constraint on the uniform \(f\) before submitting.**
   The same self-check that caught v8's bug. Plug in
   \(f = 2\cdot\mathbf 1_{[-1/4, 1/4]}\), compute \(\sigma_K(t_*)\)
   for a few \(t_*\), and verify that it *equals* the actual
   \(F(t_*) = (\text{triangle})(t_*)\) up to the Fejér-mean
   approximation error. If the constraint excludes the uniform \(f\),
   you have another bug.

3. **Extend the (W) family** from v7 to also include shifted windows
   (\(B = [\pm 1/4 - h, \pm 1/4 + h]\)) — these were originally in v4
   but got dropped. They are particularly useful for catching the
   non-symmetric extremal \(f\)'s.

If v9 does (1) and (2), it should produce the first non-trivial
*upper* bound on \(C^{\text{even}}_{6.2}\) from a Fourier-side
relaxation, and combined with v7's primal-side (W) family it should
significantly tighten the v7 numerical results.

### Side note: this v8 verdict is particularly easy to write.

I usually have to verify each step of a proof, hunt for
counterexamples, and check for inheritance. v8 hands me a complete
counterexample on a plate ("here's the admissible \(f\), here's the
inequality it violates, here's why my constraint is wrong"). The
verifier's job is just to confirm the arithmetic, which I did
above. **Thank you for the explicit counterexample.** It saved me
a lot of work, and it should be the model for future iterations.

---

**Final verdict: INVALID (self-INVALID, with the right diagnosis and the right fix planned).**
