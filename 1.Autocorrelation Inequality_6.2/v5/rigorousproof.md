# v5 — Rigorous proof of the even-restricted lower bound

## Setup

We continue with the (p, Q) lift on the f-side, but in this version we
**restrict** the relaxation to even f.  The optimum of the relaxation is
therefore a valid lower bound only on the **even-restricted constant**

$$
C^{\mathrm{even}}_{6.2} \;:=\; \inf_{f \text{ admissible AND } f \text{ even}}\, \|f * f\|_\infty.
$$

Note that $C^{\mathrm{even}}_{6.2} \ge C_{6.2}$ (restricting the f's gives
a larger inf), so a lower bound on $C^{\mathrm{even}}_{6.2}$ does **not
automatically** translate into a lower bound on $C_{6.2}$.  Nonetheless,
v5 demonstrates the relaxation technology works in principle and reproduces
the known analytic bound $C^{\mathrm{even}}_{6.2} \ge 2$ at the level of
the convex program.

## Discretization

Let $L = 1/(4N)$ and $I_j = [-1/4 + (j-1)L,\, -1/4 + jL]$ for $j = 1, \dots, 2N$.
Note that the reflection $x \mapsto -x$ exchanges $I_j$ and $I_{2N+1-j}$
because $-I_j = [(j-2N)L + 1/4 - L,\, (j-2N)L + 1/4] = I_{2N+1-j}$ (use
$L = 1/(4N)$ to verify).

Variables: $p_j = \int_{I_j} f$, $Q[j,k] = p_j p_k$, plus $\Omega \ge 0$.

## v5 constraints

In addition to v1's PSD lift, mass, marginals, $Q \ge 0$, bands, and
central window, we add:

### (E) Even-f symmetry on $p$

$$
p_j \;=\; p_{2N+1-j}, \qquad j = 1, \dots, 2N.
$$

**Proof.**  $f$ even implies $f(x) = f(-x)$, so
$$
p_{2N+1-j} \;=\; \int_{I_{2N+1-j}} f \;=\; \int_{-I_j} f
\;=\; \int_{I_j} f(-y)\,dy \;=\; \int_{I_j} f(y)\,dy \;=\; p_j. \qquad\square
$$

### (F) Even-f symmetry on $Q$

$$
Q[j, k] \;=\; Q[2N+1-j,\, 2N+1-k], \qquad j, k = 1, \dots, 2N.
$$

**Proof.**  $Q[j,k] = p_j p_k = p_{2N+1-j}\, p_{2N+1-k} = Q[2N+1-j,\, 2N+1-k]$
using (E) twice. $\square$

### (G) Anti-diagonal equals diagonal

$$
Q[j,\, 2N+1-j] \;=\; Q[j,\, j], \qquad j = 1, \dots, 2N.
$$

**Proof.**  By (E), $p_{2N+1-j} = p_j$, so
$$
Q[j,\, 2N+1-j] \;=\; p_j\, p_{2N+1-j} \;=\; p_j \cdot p_j \;=\; p_j^{2}
\;=\; Q[j, j]. \qquad\square
$$

This is the **key constraint** that breaks the diagonal-Q degeneracy of
v1–v4.  The pathological feasible point $Q = (1/(2N))\, I$ has
$Q[j, 2N+1-j] = 0$ but $Q[j,j] = 1/(2N)$, so it violates (G).  Adding (G)
forces the relaxation off this degenerate point.

## Validity of the relaxation (for even f)

Every admissible **even** f gives a feasible $(\Omega^*, p^*, Q^*)$ with
$\Omega^* = \|f * f\|_\infty$, $p_j^* = \int_{I_j} f$, and
$Q^*[j,k] = p_j^* p_k^*$:

* The v1 constraints (PSD lift, mass, marginals, $Q \ge 0$, bands, central
  window) hold by their respective derivations.
* (E), (F), (G) hold by the proofs above.

So the v5 optimum $\Omega^{(5)} \le \|f * f\|_\infty$ for every
admissible even $f$, giving
$$
\Omega^{(5)} \;\le\; C^{\mathrm{even}}_{6.2}.
$$

## Numerical results

| $N$ | $\Omega^{(5)}$ |
|-----|----------------|
| 4   | 1.800000       |
| 6   | 1.857142       |
| 8   | 1.888886       |
| 10  | 1.909089       |
| 12  | 1.923067       |

The relaxed optimum is monotonically increasing toward 2 as $N$ increases,
which matches the analytic bound $C^{\mathrm{even}}_{6.2} = 2$ achieved
by the uniform $f = 2 \cdot \mathbf 1_{[-1/4, 1/4]}$ (whose autocorrelation
is the triangle of peak 2 at $t = 0$).

## Caveat: this is a bound on the EVEN-restricted constant only

Because the constraints (E), (F), (G) all rely on $f$ being even, the v5
bound is NOT a valid lower bound on $C_{6.2}$ in general — only on
$C^{\mathrm{even}}_{6.2}$.

It is conjectured (and consistent with the SOTA $C_{6.2} \ge 1.28$) that
$C_{6.2} < 2$, so $C^{\mathrm{even}}_{6.2}$ is strictly larger than
$C_{6.2}$.  v6 will attempt to adapt the anti-diagonal-equals-diagonal
trick to a constraint that is valid for **all admissible f** by working
with both $f$ and $f \!\!\check{} (x) = f(-x)$ simultaneously and using
the fact that
$$
\int f(x) f(-x) \,dx \;=\; \frac{1}{2}\Bigl(\int f^{2} - \int (f - f \!\!\check{})^{2}/2\Bigr),
$$
which couples the evenness defect to a quantity that the relaxation can
constrain.
