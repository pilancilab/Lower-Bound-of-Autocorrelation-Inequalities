/-
  Verify/V51.lean

  Machine-checked verification of the mathematical core of the
  Autocorrelation Inequality 6.2 **v51** lower-bound convex program.

  v51 adds the "localized (Putinar) Bochner–Toeplitz" constraint using the
  localizer h(x) = cos(2πx), which is ≥ 0 on the support [-1/4, 1/4].
  For any nonnegative f supported there, ν := h·f is a positive measure, and
  its Hermitian Toeplitz matrix T_ν ⪰ 0 (Bochner). Its symbol is
      ν̂(k) = ½ ( f̂(k-1) + f̂(k+1) ).

  This file verifies the parts of that construction that a sign/index error
  would break, i.e. exactly what `convex_code.py` encodes:

    1. `sgn`, `fhat`         — the real (a,b) encoding of f̂ with sgn(0)=0.
    2. `nuhat_re`,`nuhat_im` — the real/imaginary split of ν̂ MATCHES the
                               code's `nu_re` / `nu_im` (this is the crux).
    3. `fhat_neg`,`nuhat_conj` — Hermitian symmetry ν̂(-k)=conj ν̂(k), so the
                               real block-lift `[[Re,-Im],[Im,Re]]` is valid.
    4. `cos_localizer_nonneg` — h(x)=cos(2πx) ≥ 0 on [-1/4,1/4].
    5. `localized_toeplitz_form_eq` / `_nonneg` — the localized-Bochner
       quadratic form is a nonnegative combination of squares, so EVERY
       admissible f is feasible ⇒ the constraint is a valid relaxation
       (adding it can only raise the minimum ⇒ the bound stays valid).

  The continuum Bochner theorem itself is classical; here it is proved
  directly for finite atomic measures, which is what pins down the signs.
-/
import Mathlib

open Complex BigOperators
open scoped Real

namespace V51
noncomputable section

/-! ### 1. The real (a,b) encoding of the Fourier coefficients f̂. -/

/-- Integer sign as a real number, with the convention `sgn 0 = 0`
    (this is the `_sign` helper in `convex_code.py`). -/
def sgn (m : ℤ) : ℝ := if 0 < m then 1 else if m < 0 then -1 else 0

variable (a b : ℕ → ℝ)

/-- `f̂(m) = a_{|m|} - i · sgn(m) · b_{|m|}`.  For `m ≥ 0` this is
    `a_m - i b_m`; for `m < 0` it is the conjugate `a_{|m|} + i b_{|m|}`. -/
def fhat (m : ℤ) : ℂ :=
  (a m.natAbs : ℂ) - ((sgn m * b m.natAbs : ℝ) : ℂ) * Complex.I

@[simp] lemma fhat_re (m : ℤ) : (fhat a b m).re = a m.natAbs := by
  simp [fhat]

@[simp] lemma fhat_im (m : ℤ) : (fhat a b m).im = -(sgn m * b m.natAbs) := by
  simp [fhat]

/-- The sign function flips under negation. -/
lemma sgn_neg (m : ℤ) : sgn (-m) = - sgn m := by
  rcases lt_trichotomy m 0 with h | h | h
  · have h1 : 0 < -m := by omega
    have h2 : ¬ 0 < m := by omega
    simp only [sgn, if_pos h1, if_neg h2, if_pos h]; norm_num
  · subst h; simp [sgn]
  · have h1 : ¬ 0 < -m := by omega
    have h2 : -m < 0 := by omega
    simp only [sgn, if_neg h1, if_pos h2, if_pos h]

/-- Hermitian symmetry of the underlying sequence: `f̂(-m) = conj (f̂ m)`. -/
lemma fhat_neg (m : ℤ) : fhat a b (-m) = (starRingEnd ℂ) (fhat a b m) := by
  apply Complex.ext
  · simp [fhat_re, Complex.conj_re, Int.natAbs_neg]
  · simp only [fhat_im, Complex.conj_im, Int.natAbs_neg, sgn_neg]
    ring

/-! ### 2. The localized symbol ν̂ and its real/imaginary split. -/

/-- `ν̂(k) = ½ ( f̂(k-1) + f̂(k+1) )`, the localized Bochner symbol. -/
def nuhat (k : ℤ) : ℂ := (fhat a b (k - 1) + fhat a b (k + 1)) / 2

/-- **Matches `nu_re` in `convex_code.py`:**
    `Re ν̂(k) = (a_{|k-1|} + a_{|k+1|}) / 2`. -/
theorem nuhat_re (k : ℤ) :
    (nuhat a b k).re = (a (k - 1).natAbs + a (k + 1).natAbs) / 2 := by
  simp only [nuhat, Complex.div_re, Complex.add_re, fhat_re, Complex.add_im, fhat_im,
    Complex.re_ofNat, Complex.im_ofNat, Complex.normSq_ofNat]
  ring

/-- **Matches `nu_im` in `convex_code.py`:**
    `Im ν̂(k) = -(sgn(k-1)·b_{|k-1|} + sgn(k+1)·b_{|k+1|}) / 2`. -/
theorem nuhat_im (k : ℤ) :
    (nuhat a b k).im =
      -(sgn (k - 1) * b (k - 1).natAbs + sgn (k + 1) * b (k + 1).natAbs) / 2 := by
  simp only [nuhat, Complex.div_im, Complex.add_im, fhat_im, Complex.add_re, fhat_re,
    Complex.re_ofNat, Complex.im_ofNat, Complex.normSq_ofNat]
  ring

/-- Hermitian symmetry of the Toeplitz symbol: `ν̂(-k) = conj (ν̂ k)`, so the
    matrix `T_ν[j,l] = ν̂(j-l)` is Hermitian and the real block-lift
    `[[Re, -Im], [Im, Re]] ⪰ 0` is exactly `T_ν ⪰ 0`. -/
theorem nuhat_conj (k : ℤ) :
    nuhat a b (-k) = (starRingEnd ℂ) (nuhat a b k) := by
  have h1 : (-k) - 1 = -(k + 1) := by ring
  have h2 : (-k) + 1 = -(k - 1) := by ring
  simp only [nuhat, h1, h2, fhat_neg, map_div₀, map_add, map_ofNat]
  ring

/-! ### 3. The localizer `h(x) = cos(2πx)` is nonnegative on `[-1/4, 1/4]`. -/

theorem cos_localizer_nonneg {x : ℝ} (hx : x ∈ Set.Icc (-(1/4 : ℝ)) (1/4)) :
    0 ≤ Real.cos (2 * π * x) := by
  obtain ⟨hl, hr⟩ := hx
  apply Real.cos_nonneg_of_mem_Icc
  constructor
  · nlinarith [Real.pi_pos]
  · nlinarith [Real.pi_pos]

/-! ### 4. Validity: the localized Toeplitz quadratic form is a nonnegative
        combination of squares (Bochner for finite atomic measures). -/

variable {n K : ℕ}

/-- The unit-circle Fourier atom `e^{2πix}`. -/
def atom (x : ℝ) : ℂ := Complex.exp (2 * π * Complex.I * (x : ℂ))

/-- The localized moments realized by an atomic measure `Σ_i mu_i δ_{x_i}`:
    `T[j,l] = Σ_i mu_i · conj(atom_i)^j · atom_i^l = ν̂(j-l)`.
    Here `mu_i` plays the role of `cos(2πx_i)·(weight_i) ≥ 0`. -/
def Tent (mu : Fin n → ℝ) (x : Fin n → ℝ) (j l : ℕ) : ℂ :=
  ∑ i, (mu i : ℂ) * ((starRingEnd ℂ) (atom (x i))) ^ j * (atom (x i)) ^ l

/-- The trigonometric polynomial `P_i = Σ_l w_l · atom_i^l`. -/
def Pval (x : Fin n → ℝ) (w : Fin (K + 1) → ℂ) (i : Fin n) : ℂ :=
  ∑ l : Fin (K + 1), w l * (atom (x i)) ^ (l : ℕ)

/-- **Key identity.** The Hermitian form of the localized Toeplitz matrix is a
    nonnegative combination of squared magnitudes:
    `Σ_{j,l} conj(w_j) · T[j,l] · w_l = Σ_i mu_i · |P_i|²`. -/
theorem localized_toeplitz_form_eq
    (mu : Fin n → ℝ) (x : Fin n → ℝ) (w : Fin (K + 1) → ℂ) :
    (∑ j : Fin (K + 1), ∑ l : Fin (K + 1),
        (starRingEnd ℂ) (w j) * Tent mu x j l * w l)
      = ((∑ i, mu i * Complex.normSq (Pval x w i) : ℝ) : ℂ) := by
  -- Step A: expand `Tent`, exposing the inner i-sum, giving a triple sum j,l,i.
  have stepA :
      (∑ j : Fin (K + 1), ∑ l : Fin (K + 1),
          (starRingEnd ℂ) (w j) * Tent mu x j l * w l)
        = ∑ j : Fin (K + 1), ∑ l : Fin (K + 1), ∑ i,
            (mu i : ℂ) * ((starRingEnd ℂ) (w j) * ((starRingEnd ℂ) (atom (x i))) ^ (j : ℕ))
              * (w l * (atom (x i)) ^ (l : ℕ)) := by
    refine Finset.sum_congr rfl (fun j _ => Finset.sum_congr rfl (fun l _ => ?_))
    rw [Tent, Finset.mul_sum, Finset.sum_mul]
    exact Finset.sum_congr rfl (fun i _ => by ring)
  -- Step B: swap the sums to bring i outermost:  ∑_j ∑_l ∑_i  =  ∑_i ∑_j ∑_l.
  have stepB :
      (∑ j : Fin (K + 1), ∑ l : Fin (K + 1), ∑ i,
          (mu i : ℂ) * ((starRingEnd ℂ) (w j) * ((starRingEnd ℂ) (atom (x i))) ^ (j : ℕ))
            * (w l * (atom (x i)) ^ (l : ℕ)))
        = ∑ i, ∑ j : Fin (K + 1), ∑ l : Fin (K + 1),
            (mu i : ℂ) * ((starRingEnd ℂ) (w j) * ((starRingEnd ℂ) (atom (x i))) ^ (j : ℕ))
              * (w l * (atom (x i)) ^ (l : ℕ)) := by
    have inner : ∀ j : Fin (K + 1),
        (∑ l : Fin (K + 1), ∑ i,
            (mu i : ℂ) * ((starRingEnd ℂ) (w j) * ((starRingEnd ℂ) (atom (x i))) ^ (j : ℕ))
              * (w l * (atom (x i)) ^ (l : ℕ)))
          = ∑ i, ∑ l : Fin (K + 1),
            (mu i : ℂ) * ((starRingEnd ℂ) (w j) * ((starRingEnd ℂ) (atom (x i))) ^ (j : ℕ))
              * (w l * (atom (x i)) ^ (l : ℕ)) := fun j => Finset.sum_comm
    rw [Finset.sum_congr rfl (fun j _ => inner j), Finset.sum_comm]
  -- Step C: each i-term factors as `mu_i · conj(P_i) · P_i = mu_i · |P_i|²`.
  have stepC : ∀ i : Fin n,
      (∑ j : Fin (K + 1), ∑ l : Fin (K + 1),
          (mu i : ℂ) * ((starRingEnd ℂ) (w j) * ((starRingEnd ℂ) (atom (x i))) ^ (j : ℕ))
            * (w l * (atom (x i)) ^ (l : ℕ)))
        = ((mu i * Complex.normSq (Pval x w i) : ℝ) : ℂ) := by
    intro i
    have hconjP : (starRingEnd ℂ) (Pval x w i)
        = ∑ j : Fin (K + 1), (starRingEnd ℂ) (w j) * ((starRingEnd ℂ) (atom (x i))) ^ (j : ℕ) := by
      rw [Pval, map_sum]
      exact Finset.sum_congr rfl (fun j _ => by rw [map_mul, map_pow])
    calc
      (∑ j : Fin (K + 1), ∑ l : Fin (K + 1),
          (mu i : ℂ) * ((starRingEnd ℂ) (w j) * ((starRingEnd ℂ) (atom (x i))) ^ (j : ℕ))
            * (w l * (atom (x i)) ^ (l : ℕ)))
          = (mu i : ℂ) *
              ((∑ j : Fin (K + 1), (starRingEnd ℂ) (w j) * ((starRingEnd ℂ) (atom (x i))) ^ (j : ℕ))
                * (∑ l : Fin (K + 1), w l * (atom (x i)) ^ (l : ℕ))) := by
            rw [Finset.sum_mul_sum, Finset.mul_sum]
            refine Finset.sum_congr rfl (fun j _ => ?_)
            rw [Finset.mul_sum]
            exact Finset.sum_congr rfl (fun l _ => by ring)
      _ = (mu i : ℂ) * ((starRingEnd ℂ) (Pval x w i) * Pval x w i) := by
            rw [← hconjP]; rfl
      _ = ((mu i * Complex.normSq (Pval x w i) : ℝ) : ℂ) := by
            rw [mul_comm ((starRingEnd ℂ) (Pval x w i)) (Pval x w i), Complex.mul_conj]
            push_cast; ring
  rw [stepA, stepB, Complex.ofReal_sum]
  exact Finset.sum_congr rfl (fun i _ => stepC i)

/-- **Validity of the v51 constraint.** For nonnegative weights `mu_i ≥ 0`
    (in the application `mu_i = cos(2πx_i)·f_i ≥ 0` on the support), the
    localized Toeplitz Hermitian form is real and `≥ 0` — i.e. `T_ν ⪰ 0`
    holds at every genuine (atomic) `f`, so the constraint cuts off no
    admissible point. -/
theorem localized_toeplitz_form_nonneg
    (mu : Fin n → ℝ) (x : Fin n → ℝ) (hmu : ∀ i, 0 ≤ mu i)
    (w : Fin (K + 1) → ℂ) :
    0 ≤ (∑ j : Fin (K + 1), ∑ l : Fin (K + 1),
          (starRingEnd ℂ) (w j) * Tent mu x j l * w l).re := by
  rw [localized_toeplitz_form_eq, Complex.ofReal_re]
  apply Finset.sum_nonneg
  intro i _
  exact mul_nonneg (hmu i) (Complex.normSq_nonneg _)

end
end V51
