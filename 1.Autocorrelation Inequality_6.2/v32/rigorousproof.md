# v32 — Rigorous proof: v30 + v11 windows

## Setup

v32 combines v30's pointwise Fejér mean tests with v11's shifted
windowed (W) constraints on the (p, Q) lift.

## Validity

Each constraint is individually valid for any admissible $f$:
* v14 autocorrelation Fejér (proven in v14)
* v30 pointwise Fejér mean tests at multiple $t_*$ (proven in v30)
* v11 shifted (W) windowed integral inequalities (proven in v7/v11)
* (p, Q) PSD lift, mass, marginals, $Q \ge 0$ (proven in v0–v7)

The conjunction is valid; the rank-one truth from any admissible $f$
satisfies all of them simultaneously.

## Numerical observation

At $N = 16$, $K = 16$, $n_\text{tests} = 33$:
$\Omega^{(32)} = 1.147$, **identical** to v30 at the same parameters.

The (W) windowed constraints **do not improve** the v30 bound at
this scale — they are dominated by v14's autocorrelation Fejér +
v30's pointwise Fejér.

## Conclusion

v32 is rigorously valid but produces no numerical improvement over
v30.  The dominant constraint at the LP optimum is the v30
pointwise Fejér mean test.

The best rigorous bound on $C_{6.2}$ remains v30/v31's
$\Omega \approx 1.267$ at $(N, K) = (2000, 32)$, $0.013$ below
the SOTA $C_{6.2} \ge 1.28$.

## Project summary

The project has produced rigorous lower bounds via convex relaxation
across 30+ iterations.  The major milestones are:

* **v7 (VALID, even-restricted)**: First clean rigorous reboot.
  $\Omega \approx 1.34$.
* **v9 (VALID, even-restricted)**: Fejér kernel test functional.
  $\Omega \approx 1.69$ at large $N, K$.
* **v14 (VALID, unrestricted)**: Autocorrelation Fejér via
  $|\hat F(k)|^2$ lift.  First non-trivial bound on the
  unrestricted $C_{6.2}$.
* **v30 / v31 (VALID, unrestricted)**: v14 + pointwise Fejér mean
  tests at multiple $t_*$.  $\Omega \approx 1.267$ at $N = 2000$,
  $K = 32$.

The unrestricted bound has reached within $\sim 0.01$ of the SOTA.
Further work may push past it via:
1. Larger numerical computation ($N = 10^4, K = 64$).
2. Additional Fourier-side constraints (Bochner Toeplitz) at higher
   moments.
3. Higher-order Lasserre lifts on the moment matrix $M$.

Each direction is rigorous in principle, but requires either compute
or new mathematical insights.
