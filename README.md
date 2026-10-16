# Lower-Bound-of-Autocorrelation-Inequalities

Research workspace for developing and checking optimization bounds for convolution inequalities and the Erdős minimum overlap problem. It combines versioned mathematical proposals, Python optimization programs, numerical verification artifacts, and Lean formalizations.

The main workflow is to propose a constraint, justify it mathematically, review the proof, implement the relaxation, and evaluate the result. A separate constructive search produces feasible examples for upper bounds.

## Project layout

| Directory | Contents |
| --- | --- |
| [1.Autocorrelation Inequality_6.2](<1.Autocorrelation Inequality_6.2/>) | Lower-bound experiments for the self-convolution inequality, with versions `v0`–`v54`, a `v51_verified` variant using analytic cell extrema, and objective plots. |
| [2.Erdos minimum overlap problem](<2.Erdos minimum overlap problem/>) | Convex lower-bound experiments `v0`–`v7`, based on Fourier, moment, and positivity constraints. |
| [3.Upper bound trial](<3.Upper bound trial/>) | Constructive upper-bound searches `v0`–`v15`, saved step functions, independent numerical verifiers, and reference constructions. |
| [Convex_Problems_Verified](Convex_Problems_Verified/) | Shared convex models, dual-witness certification routines, regional verification scripts, branch-and-bound experiments, and saved results. |
| [Lean_Verification](Lean_Verification/README.md) | Lean 4/Mathlib formalizations of selected constraint identities and soundness arguments. |
| [White-Minimum-Overlap](White-Minimum-Overlap/README.md) | Reference OPL programs, feasible solutions, and a local copy of the paper motivating the minimum-overlap relaxations. |

This is a collection of research scripts and artifacts. There is no root-level package, unified command-line interface, or pinned Python dependency manifest.

## Mathematical problems

**Autocorrelation Inequality 6.2.** Despite the directory name, the quantity studied is the self-convolution `(f * f)(t)`. The goal is a lower bound on the largest constant `C_6.2` such that every admissible nonnegative function satisfies

```text
max_{|t| ≤ 1/2} (f * f)(t) ≥ C_6.2 (∫_{-1/4}^{1/4} f(x) dx)².
```

See the [problem statement](<1.Autocorrelation Inequality_6.2/Problem.md>).

**Erdős minimum overlap.** The implemented searches use functions `f, g : [-1,1] → [0,1]` with `g = 1 - f` on that interval and `∫ f = 1`. Both functions are extended by zero outside the interval. The objective is

```text
inf_f sup_{x ∈ [-2,2]} ∫_{-1}^{1} f(t) g(t+x) dt.
```

Convex relaxations seek lower bounds. Explicit feasible step functions give upper bounds. For a vector `a` on `n` equal bins, feasibility requires `0 ≤ a_i ≤ 1` and `sum(a) = n/2`; overlap values are computed by linear cross-correlation with `1-a`, scaled by `2/n`.

See the [lower-bound problem notes](<2.Erdos minimum overlap problem/Problem.md>) and [constructive-search specification](<3.Upper bound trial/proposal.md>). The mass normalization above is essential to the nontrivial minimum-overlap problem.

## Python setup

Use Python 3.10 or newer for the syntax used across the scripts. From the workspace root:

```sh
python3 -m venv .venv
. .venv/bin/activate
python -m pip install numpy scipy cvxpy scs mpmath matplotlib
```

Most SDP experiments default to MOSEK. To use those paths, install its Python package and configure a valid MOSEK license:

```sh
python -m pip install Mosek
```

Some experiment scripts support SCS; the dual-witness solver in `Convex_Problems_Verified/certified_bound.py` explicitly uses MOSEK. PDF table generation additionally imports `fpdf`. The reference OPL programs require IBM CPLEX/OPL separately.

## Running examples

Run these commands from the workspace root with the Python environment activated. Quote paths containing spaces.

**Check a saved upper-bound construction.** This reads the existing sequence, checks box and mass feasibility, and recomputes overlap values using FFT and, for this sequence size, a direct sum:

```sh
python "3.Upper bound trial/v15/verifier.py"
```

To check another saved sequence with the same verifier:

```sh
python "3.Upper bound trial/v15/verifier.py" "3.Upper bound trial/v14/sequence.json"
```

**Run a small autocorrelation relaxation.** Use the analytic-cell-bound variant with modest parameters before attempting larger solves:

```sh
(
  cd "1.Autocorrelation Inequality_6.2/v51_verified"
  python - <<'PY'
from convex_code import AutocorrLowerBoundV51v

problem = AutocorrLowerBoundV51v(N=8, K=8)
print(problem.solve(solver="MOSEK"))
PY
)
```

**Run dual certification at a small size.** This writes a result JSON and a dual witness under `Convex_Problems_Verified/y_witnesses/`:

```sh
python Convex_Problems_Verified/verify_prob62.py \
  --N 8 --K 8 \
  --out Convex_Problems_Verified/prob62_N8_K8_results.json
```

The defaults are `N=200, K=96`; the saved run at those settings took approximately 41 minutes. Inspect `valid` and any `error` in the result JSON as well as the printed output.

**Run a constructive search.** Version 15 starts from the saved `v13` sequence and searches for improvements with exact integer pair moves:

```sh
BUDGET_SEC=60 python "3.Upper bound trial/v15/solver.py"
```

This overwrites `3.Upper bound trial/v15/sequence.json`. The budget controls the search phase; setup and final evaluation add runtime. Copy the version directory first if you want to preserve its saved result.

Larger regional searches live in `Convex_Problems_Verified/bnb/`. Their grid sizes, worker counts, thresholds, and output paths are configured in the individual scripts.

## Lean verification

The Lean project pins `leanprover/lean4:v4.32.0` and Mathlib `v4.32.0`. With Elan/Lake installed:

```sh
(
  export PATH="$HOME/.elan/bin:$PATH"
  cd Lean_Verification
  lake exe cache get
  lake build
)
```

See the [Lean README](Lean_Verification/README.md) for theorem coverage and limitations, and [convex_programs.tex](Lean_Verification/convex_programs.tex) for the transcribed optimization programs. The formalizations cover selected identities and soundness arguments, including atomic-measure models; they do not certify the numerical SDP solves. Some analytic ingredients remain cited rather than formalized.

## Reading results and extending experiments

Most lower-bound version directories contain `convex_code.py`, `rigorousproof.md`, and `verdict.md`. Read the proof and verdict together before reusing a proposed constraint. A later version number alone does not establish a better bound.

Upper-bound directories contain `solver.py`, `verifier.py`, and `sequence.json`. Shared certification outputs include result JSON files, CSV tables, plots, and compressed `.npz` dual witnesses.

Two examples already stored in the workspace are:

| Artifact | Recorded result |
| --- | --- |
| [Autocorrelation certification JSON](Convex_Problems_Verified/prob62_N200_K96_results.json) | `N=200`, `K=96`, reported lower bound `1.293736833333943`, `valid: true`. |
| [Version 15 upper-bound sequence](<3.Upper bound trial/v15/sequence.json>) | `n=600`, reported upper bound `0.38087031046756176`, plus an exact rational objective string. |

These are recorded outputs, not claims about current literature records. Distinguish a proof-review verdict, a solver's numerical feasibility check, a dual certificate, and a Lean theorem when reporting results. The upper-bound verifier uses floating-point tolerances. Dual certification depends on the supplied coefficient-error bounds and an a priori primal norm bound; `verify_prob62.py` supplies coefficient-error budgets of `1e-12`.

For new experiments, preserve the previous version and its artifacts, record the mathematical change and its justification, review the constraint, then save the code, parameters, solver status, feasibility information, and resulting bounds together.
