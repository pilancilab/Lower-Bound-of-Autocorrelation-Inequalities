# Lean verification of two lower-bound convex programs

Machine-checked (Lean 4 + Mathlib, toolchain `leanprover/lean4:v4.32.0`) core of:

1. **v51** — Autocorrelation Inequality 6.2, localized (Putinar) Bochner
   constraint with localizer `h(x)=cos(2πx)`.  → `Verify/V51.lean`
2. **Erdős minimum overlap, White v7** — the Dirichlet-kernel Bochner
   "sandwich" constraint.  → `Verify/ErdosV7.lean`

The exact convex programs (transcribed from the `convex_code.py` files) are in
[`convex_programs.tex`](convex_programs.tex).

## Build / reproduce

```sh
export PATH="$HOME/.elan/bin:$PATH"
cd Lean_Verification
lake exe cache get      # prebuilt Mathlib oleans (first time only)
lake build              # compiles Verify.V51 and Verify.ErdosV7
```

Both files compile with **no `sorry`/`admit`**; every headline theorem depends
only on the three standard axioms `propext`, `Classical.choice`, `Quot.sound`
(check with `#print axioms <name>`).

## What is proved

### v51 (`Verify/V51.lean`)
- `nuhat_re`, `nuhat_im` — the real/imag split of the localized symbol
  `ν̂(k)=½(f̂(k-1)+f̂(k+1))` **matches** the code's `nu_re`/`nu_im`
  (`sgn(0)=0`, indices `|k±1|`).  *This is the crux sign/index check.*
- `fhat_neg`, `nuhat_conj` — Hermitian symmetry, justifying the real block lift.
- `cos_localizer_nonneg` — `cos(2πx) ≥ 0` on `[-1/4,1/4]`.
- `localized_toeplitz_form_eq/_nonneg` — the localized Toeplitz Hermitian form
  is a nonnegative sum of squared magnitudes ⇒ the constraint is a **valid
  relaxation** (cuts off no admissible `f`).

### v51 FULL program soundness (`Verify/V51Program.lean`)
Models an admissible `f` as a nonnegative atomic measure and checks **every**
constraint block of `prob6.2.py` holds at the truth (⇒ valid lower bound):

| Block | meaning | status |
|---|---|---|
| P1 mass, `a₀=1`, `b₀=0` | normalization | ✓ `mass_eq_one`, `b_zero` |
| P2 `M ⪰ 0` + first row | moment lift | ✓ `outer_nonneg`, `outer_quad` |
| P3 cell bounds | `cos/sin` cell extrema | ✓ `cos_min_attained`, `sin_min_attained` |
| P6 plain Bochner `T_f ⪰ 0` | `f ≥ 0` | ✓ `plain_bochner_form_nonneg` |
| P7 localized Bochner | v51 new | ✓ (in `V51.lean`) |
| P4 `v14` energy | `Ω ≥ 1+2Σwₖvₖ` | encoding ✓ `v14_diag_sum`; bound = cited |
| P5 Fejér–Riesz | `Ω ≥ supₜ σ_K` | encoding ✓ `fejer_coeff_re/im`; bound = cited |

P4/P5's `Ω`-lower-bound direction uses two classical analytic theorems
(**Fejér–Riesz spectral factorization** and the **Fejér-kernel peak bound**
`σ_K^{f*f} ≤ ‖f*f‖_∞`). These are named in the file and left cited — the code
only needs to *encode* them, which the coefficient identities certify.
Numerically, the P3 analytic extrema were also stress-tested vs a fine grid
(30k random cases, 0 unsafe cases).

### Erdős v7 (`Verify/ErdosV7.lean`)
- `dirichlet_cos_integral`, `dirichlet_sin_integral` — the indicator integral
  `∫_{-1}^1 e^{-iπkx/2}dx = 4·D(k)` is real, `D(k)=sin(kπ/2)/(kπ)`
  (source of the constant matrix `S`).
- `sandwich_g0`, `sandwich_fext` — the entrywise identities
  `T^{g0}=½I+S-T^f` and `T^{fext}=½I-S+T^f` (**exact signs** of the two LMIs).
- `g0_nonneg`, `fext_nonneg` — the completions are `≥ 0` on `[-2,2]`.
- `bochner_toeplitz_form_nonneg` — Bochner for atomic measures (PSD mechanism).

## Scope

The proofs certify **soundness** of each *new* constraint (every admissible `f`
stays feasible ⇒ the program's optimum is a valid lower bound). They do **not**
re-run or certify the numerical SDP solve. No sign/index/factor error was found
in either new constraint.
