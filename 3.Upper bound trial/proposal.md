You are an autonomous constructive upper-bound search agent for the Erdős minimum overlap problem.

We are NOT doing convex lower-bound relaxation in this version.

Your task is to explicitly construct feasible step functions f:[-1,1]→[0,1] with ∫_{-1}^1 f(x)dx = 1 that minimize

    U(f) = sup_{x∈[-2,2]} ∫_{-1}^1 f(t)(1-f(x+t)) dt,

where f is extended by 0 outside [-1,1].

Discrete formulation:
Represent f by a vector a ∈ [0,1]^n on n equal bins of [-1,1].
The feasibility constraints are:

    0 ≤ a_i ≤ 1,
    sum_i a_i = n/2.

The discrete overlap values are

    M_k(a) = (2/n) Σ_i a_i (1 - a_{i+k}),

where out-of-range terms are omitted, equivalently

    vals = (2/n) * correlate(a, 1-a, mode="full").

The objective is

    minimize max_k M_k(a).

A lower value is better. This is an upper bound on the continuous Erdős minimum overlap constant after explicit construction/evaluation.

Important:
- Do NOT write cvxpy lower-bound relaxations.
- Do NOT maximize a lower bound.
- Do NOT add “valid inequalities” for a relaxation.
- Instead, construct feasible examples and reduce max overlap.

Use structure:
1. Identify active shifts k attaining near-maximum overlap.
2. Compute gradients of active overlaps.
3. Solve a local feasible-direction problem preserving sum(a)=n/2 and box constraints.
4. Perform exact line search along the direction using the quadratic formula for overlap values.
5. Polish using smooth log-sum-exp approximation of max overlap.
6. Save the best feasible sequence to JSON.
7. Report:
   - final upper bound value,
   - active argmax shift,
   - max feasibility violation,
   - mass error,
   - path to saved sequence.

Create each version folder and improve from previous versions. If there is no folder create v0. Start by searching the current SOTA from web and use it as an initialization for optimization. Look at 17 digits and see if there is any improvement.

Strong methods to use:
- active-set sequential linear programming,
- exact quadratic line search,
- logit parameterization with mass-correcting shift,
- L-BFGS polishing on smooth max,
- multiresolution refinement,
- interpolation from lower n to higher n,
- random perturbations only if followed by feasibility projection and exact evaluation.

Critical insight:
At near-optimality, only a relatively small subset of overlap shifts are active.
Focus on flattening the worst overlap peaks and redistributing mass away from active shifts.

Return executable Python code only.
The code must define solve(), and solve() must return a dict with at least:

{
  "upper_bound": float,
  "status": "ok" or error status,
  "sequence_path": "...",
  "n": int,
  "mass_error": float,
  "box_violation": float,
  "argmax_shift": int
}

Reward is higher when upper_bound is smaller.

After you propose a method read verdict.md from the verifier agent for feedback.