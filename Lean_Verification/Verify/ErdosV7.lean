/-
  Verify/ErdosV7.lean

  Machine-checked verification of the mathematical core of the **v7**
  lower-bound convex program for the Erdős minimum-overlap problem
  ("Dirichlet-kernel Bochner sandwich", White's period-4 setup).

  On the period-4 torus [-2,2], f is supported in [-1,1] with 0 ≤ f ≤ 1 and
  ∫_{-1}^1 f = 1 (so the 0-th Fourier coefficient is 1/4).  v7 adds two LMIs
  from Bochner–Herglotz applied to two nonnegative completions of f:

      g_0     := (1-f)·1_{[-1,1]}                         (≥ 0 on [-2,2])
      f_ext   := f·1_{[-1,1]} + 1_{[-2,2]} - 1_{[-1,1]}   (≥ 0 on [-2,2])

  With the Dirichlet-kernel constant matrix S_{kl}=D(k-l),
  D(j)=sin(jπ/2)/(jπ):

      T^{g_0}   = ½ I + S - T^f      (⪰ 0)
      T^{f_ext} = ½ I - S + T^f      (⪰ 0)

  Verified here (the parts a sign/factor/index error would break):
    1. `dirichlet_cos_integral` / `dirichlet_sin_integral`
         `∫_{-1}^1 cos(πkx/2)dx = 4·D(k)` and `∫_{-1}^1 sin(πkx/2)dx = 0`,
         i.e. the indicator integral `∫_{-1}^1 e^{-iπkx/2}dx = 4 D(k)` is real
         and equals the `S_mat` entries in `convex_code.py`.
    2. `sandwich_g0` / `sandwich_fext`
         the entrywise Toeplitz identities `T^{g_0}=½I+S-T^f`,
         `T^{f_ext}=½I-S+T^f` (fixing the signs of `half_I_S ± S_block ∓ Sf`).
    3. `g0_nonneg` / `fext_nonneg` — the completions are ≥ 0 on [-2,2].
    4. `bochner_toeplitz_form_nonneg` — a nonnegative (atomic) measure yields a
         PSD Toeplitz form: the mechanism giving `T^{g_0},T^{f_ext} ⪰ 0`.
-/
import Mathlib

open Complex intervalIntegral BigOperators
open scoped Real

namespace ErdosV7
noncomputable section

/-! ### 1. The Dirichlet-kernel constant and the indicator integral. -/

/-- Dirichlet constant `D(j) = sin(jπ/2)/(jπ)`, with `D 0 = 0`. -/
def D (j : ℤ) : ℝ := if j = 0 then 0 else Real.sin (j * π / 2) / (j * π)

@[simp] lemma D_zero : D 0 = 0 := by simp [D]

/-- Antiderivative computation `∫_{-1}^1 cos(c x) dx = 2 sin c / c`. -/
private theorem integral_cos_lin (c : ℝ) (hc : c ≠ 0) :
    (∫ x in (-1 : ℝ)..1, Real.cos (c * x)) = 2 * Real.sin c / c := by
  have hderiv : ∀ x ∈ Set.uIcc (-1 : ℝ) 1,
      HasDerivAt (fun x => Real.sin (c * x) / c) (Real.cos (c * x)) x := by
    intro x _
    have hg : HasDerivAt (fun x : ℝ => c * x) c x := by
      simpa using (hasDerivAt_id x).const_mul c
    have hs : HasDerivAt (fun x : ℝ => Real.sin (c * x)) (Real.cos (c * x) * c) x :=
      (Real.hasDerivAt_sin (c * x)).comp x hg
    have hd := hs.div_const c
    rwa [mul_div_cancel_right₀ _ hc] at hd
  have hcont : Continuous (fun x : ℝ => Real.cos (c * x)) := by fun_prop
  rw [integral_eq_sub_of_hasDerivAt hderiv (hcont.intervalIntegrable _ _)]
  rw [mul_one, mul_neg_one, Real.sin_neg]
  ring

/-- Antiderivative computation `∫_{-1}^1 sin(c x) dx = 0` (odd integrand). -/
private theorem integral_sin_lin (c : ℝ) (hc : c ≠ 0) :
    (∫ x in (-1 : ℝ)..1, Real.sin (c * x)) = 0 := by
  have hderiv : ∀ x ∈ Set.uIcc (-1 : ℝ) 1,
      HasDerivAt (fun x => -Real.cos (c * x) / c) (Real.sin (c * x)) x := by
    intro x _
    have hg : HasDerivAt (fun x : ℝ => c * x) c x := by
      simpa using (hasDerivAt_id x).const_mul c
    have hcos : HasDerivAt (fun x : ℝ => Real.cos (c * x)) (-Real.sin (c * x) * c) x :=
      (Real.hasDerivAt_cos (c * x)).comp x hg
    have hd := (hcos.neg).div_const c
    rw [neg_mul, neg_neg, mul_div_cancel_right₀ _ hc] at hd
    simpa [neg_div] using hd
  have hcont : Continuous (fun x : ℝ => Real.sin (c * x)) := by fun_prop
  rw [integral_eq_sub_of_hasDerivAt hderiv (hcont.intervalIntegrable _ _)]
  rw [mul_one, mul_neg_one, Real.cos_neg]
  ring

set_option maxHeartbeats 1000000 in
/-- **The real part of the indicator integral** equals `4·D(k)`:
    `∫_{-1}^1 cos(πkx/2) dx = 4·sin(kπ/2)/(kπ)` for `k ≠ 0`. -/
theorem dirichlet_cos_integral (k : ℤ) (hk : k ≠ 0) :
    (∫ x in (-1 : ℝ)..1, Real.cos (π * (k : ℝ) * x / 2))
      = 4 * Real.sin ((k : ℝ) * π / 2) / ((k : ℝ) * π) := by
  have hkR : (k : ℝ) ≠ 0 := by exact_mod_cast hk
  have hc : π * (k : ℝ) / 2 ≠ 0 := by
    apply div_ne_zero (mul_ne_zero Real.pi_ne_zero hkR); norm_num
  have hrw : (∫ x in (-1 : ℝ)..1, Real.cos (π * (k : ℝ) * x / 2))
      = ∫ x in (-1 : ℝ)..1, Real.cos ((π * (k : ℝ) / 2) * x) :=
    intervalIntegral.integral_congr (fun x _ => by congr 1; ring)
  rw [hrw, integral_cos_lin (π * (k : ℝ) / 2) hc]
  have hs : Real.sin (π * (k : ℝ) / 2) = Real.sin ((k : ℝ) * π / 2) := by congr 1; ring
  rw [hs, div_div_eq_mul_div]
  ring

/-- The imaginary part of the indicator integral vanishes:
    `∫_{-1}^1 sin(πkx/2) dx = 0`, so the Dirichlet constant `D` is real. -/
theorem dirichlet_sin_integral (k : ℤ) (hk : k ≠ 0) :
    (∫ x in (-1 : ℝ)..1, Real.sin (π * (k : ℝ) * x / 2)) = 0 := by
  have hkR : (k : ℝ) ≠ 0 := by exact_mod_cast hk
  have hc : π * (k : ℝ) / 2 ≠ 0 := by
    apply div_ne_zero (mul_ne_zero Real.pi_ne_zero hkR); norm_num
  have hrw : (∫ x in (-1 : ℝ)..1, Real.sin (π * (k : ℝ) * x / 2))
      = ∫ x in (-1 : ℝ)..1, Real.sin ((π * (k : ℝ) / 2) * x) :=
    intervalIntegral.integral_congr (fun x _ => by congr 1; ring)
  rw [hrw, integral_sin_lin (π * (k : ℝ) / 2) hc]

/-! ### 2. The Toeplitz sandwich identities (sign check on `½I ± S ∓ T^f`). -/

variable (fhat : ℤ → ℂ)

/-- Fourier coefficients of `g_0 = (1-f)·1_{[-1,1]}`:
    `ĝ_0(0)=1/4`, `ĝ_0(d)=D(d)-f̂(d)` for `d≠0`. -/
def ghat0 (d : ℤ) : ℂ := if d = 0 then (1 / 4 : ℂ) else (D d : ℂ) - fhat d

/-- Fourier coefficients of `f_ext = f·1_{[-1,1]} + 1_{[-2,2]} - 1_{[-1,1]}`:
    `f̂_ext(0)=3/4`, `f̂_ext(d)=f̂(d)-D(d)` for `d≠0`. -/
def fexthat (d : ℤ) : ℂ := if d = 0 then (3 / 4 : ℂ) else fhat d - (D d : ℂ)

/-- **Sandwich identity for `g_0`:** the Toeplitz symbol of `g_0` equals the
    entry of `½ I + S - T^f`. (Diagonal: `1/4 = 1/2 + 0 - 1/4`.) -/
theorem sandwich_g0 (h0 : fhat 0 = 1 / 4) (k l : ℤ) :
    ghat0 fhat (k - l)
      = (if k = l then (1 / 2 : ℂ) else 0) + (D (k - l) : ℂ) - fhat (k - l) := by
  unfold ghat0
  by_cases h : k = l
  · have hkl : k - l = 0 := sub_eq_zero_of_eq h
    rw [hkl, if_pos h]
    simp only [if_pos (rfl : (0 : ℤ) = 0), D_zero, Complex.ofReal_zero, h0]
    norm_num
  · have hd : k - l ≠ 0 := sub_ne_zero.mpr h
    rw [if_neg hd, if_neg h]
    ring

/-- **Sandwich identity for `f_ext`:** the Toeplitz symbol of `f_ext` equals the
    entry of `½ I - S + T^f`. (Diagonal: `3/4 = 1/2 - 0 + 1/4`.) -/
theorem sandwich_fext (h0 : fhat 0 = 1 / 4) (k l : ℤ) :
    fexthat fhat (k - l)
      = (if k = l then (1 / 2 : ℂ) else 0) - (D (k - l) : ℂ) + fhat (k - l) := by
  unfold fexthat
  by_cases h : k = l
  · have hkl : k - l = 0 := sub_eq_zero_of_eq h
    rw [hkl, if_pos h]
    simp only [if_pos (rfl : (0 : ℤ) = 0), D_zero, Complex.ofReal_zero, h0]
    norm_num
  · have hd : k - l ≠ 0 := sub_ne_zero.mpr h
    rw [if_neg hd, if_neg h]
    ring

/-! ### 3. The completions are nonnegative on `[-2,2]` (the Bochner premise). -/

variable (f : ℝ → ℝ)

/-- `g_0 = (1-f)` on `[-1,1]`, `0` outside. -/
def g0 (x : ℝ) : ℝ := if |x| ≤ 1 then 1 - f x else 0

/-- `f_ext = f` on `[-1,1]`, `1` outside. -/
def fext (x : ℝ) : ℝ := if |x| ≤ 1 then f x else 1

theorem g0_nonneg (h1 : ∀ x, |x| ≤ 1 → f x ≤ 1) (x : ℝ) : 0 ≤ g0 f x := by
  unfold g0
  split_ifs with h
  · linarith [h1 x h]
  · norm_num

theorem fext_nonneg (h0 : ∀ x, |x| ≤ 1 → 0 ≤ f x) (x : ℝ) : 0 ≤ fext f x := by
  unfold fext
  split_ifs with h
  · exact h0 x h
  · norm_num

/-! ### 4. Bochner for finite atomic measures: PSD of the Toeplitz form. -/

variable {n K : ℕ}

/-- Toeplitz matrix of the atomic measure `Σ_i mu_i δ` with unit atoms `z_i`:
    `T[j,l] = Σ_i mu_i · conj(z_i)^j · z_i^l`. -/
def ToeplitzZ (mu : Fin n → ℝ) (z : Fin n → ℂ) (j l : ℕ) : ℂ :=
  ∑ i, (mu i : ℂ) * ((starRingEnd ℂ) (z i)) ^ j * (z i) ^ l

/-- Trigonometric polynomial `P_i = Σ_l w_l · z_i^l`. -/
def PZ (z : Fin n → ℂ) (w : Fin (K + 1) → ℂ) (i : Fin n) : ℂ :=
  ∑ l : Fin (K + 1), w l * (z i) ^ (l : ℕ)

/-- The Hermitian form is a nonnegative combination of squared magnitudes. -/
theorem bochner_toeplitz_form_eq
    (mu : Fin n → ℝ) (z : Fin n → ℂ) (w : Fin (K + 1) → ℂ) :
    (∑ j : Fin (K + 1), ∑ l : Fin (K + 1),
        (starRingEnd ℂ) (w j) * ToeplitzZ mu z j l * w l)
      = ((∑ i, mu i * Complex.normSq (PZ z w i) : ℝ) : ℂ) := by
  have stepA :
      (∑ j : Fin (K + 1), ∑ l : Fin (K + 1),
          (starRingEnd ℂ) (w j) * ToeplitzZ mu z j l * w l)
        = ∑ j : Fin (K + 1), ∑ l : Fin (K + 1), ∑ i,
            (mu i : ℂ) * ((starRingEnd ℂ) (w j) * ((starRingEnd ℂ) (z i)) ^ (j : ℕ))
              * (w l * (z i) ^ (l : ℕ)) := by
    refine Finset.sum_congr rfl (fun j _ => Finset.sum_congr rfl (fun l _ => ?_))
    rw [ToeplitzZ, Finset.mul_sum, Finset.sum_mul]
    exact Finset.sum_congr rfl (fun i _ => by ring)
  have stepB :
      (∑ j : Fin (K + 1), ∑ l : Fin (K + 1), ∑ i,
          (mu i : ℂ) * ((starRingEnd ℂ) (w j) * ((starRingEnd ℂ) (z i)) ^ (j : ℕ))
            * (w l * (z i) ^ (l : ℕ)))
        = ∑ i, ∑ j : Fin (K + 1), ∑ l : Fin (K + 1),
            (mu i : ℂ) * ((starRingEnd ℂ) (w j) * ((starRingEnd ℂ) (z i)) ^ (j : ℕ))
              * (w l * (z i) ^ (l : ℕ)) := by
    have inner : ∀ j : Fin (K + 1),
        (∑ l : Fin (K + 1), ∑ i,
            (mu i : ℂ) * ((starRingEnd ℂ) (w j) * ((starRingEnd ℂ) (z i)) ^ (j : ℕ))
              * (w l * (z i) ^ (l : ℕ)))
          = ∑ i, ∑ l : Fin (K + 1),
            (mu i : ℂ) * ((starRingEnd ℂ) (w j) * ((starRingEnd ℂ) (z i)) ^ (j : ℕ))
              * (w l * (z i) ^ (l : ℕ)) := fun j => Finset.sum_comm
    rw [Finset.sum_congr rfl (fun j _ => inner j), Finset.sum_comm]
  have stepC : ∀ i : Fin n,
      (∑ j : Fin (K + 1), ∑ l : Fin (K + 1),
          (mu i : ℂ) * ((starRingEnd ℂ) (w j) * ((starRingEnd ℂ) (z i)) ^ (j : ℕ))
            * (w l * (z i) ^ (l : ℕ)))
        = ((mu i * Complex.normSq (PZ z w i) : ℝ) : ℂ) := by
    intro i
    have hconjP : (starRingEnd ℂ) (PZ z w i)
        = ∑ j : Fin (K + 1), (starRingEnd ℂ) (w j) * ((starRingEnd ℂ) (z i)) ^ (j : ℕ) := by
      rw [PZ, map_sum]
      exact Finset.sum_congr rfl (fun j _ => by rw [map_mul, map_pow])
    calc
      (∑ j : Fin (K + 1), ∑ l : Fin (K + 1),
          (mu i : ℂ) * ((starRingEnd ℂ) (w j) * ((starRingEnd ℂ) (z i)) ^ (j : ℕ))
            * (w l * (z i) ^ (l : ℕ)))
          = (mu i : ℂ) *
              ((∑ j : Fin (K + 1), (starRingEnd ℂ) (w j) * ((starRingEnd ℂ) (z i)) ^ (j : ℕ))
                * (∑ l : Fin (K + 1), w l * (z i) ^ (l : ℕ))) := by
            rw [Finset.sum_mul_sum, Finset.mul_sum]
            refine Finset.sum_congr rfl (fun j _ => ?_)
            rw [Finset.mul_sum]
            exact Finset.sum_congr rfl (fun l _ => by ring)
      _ = (mu i : ℂ) * ((starRingEnd ℂ) (PZ z w i) * PZ z w i) := by
            rw [← hconjP]; rfl
      _ = ((mu i * Complex.normSq (PZ z w i) : ℝ) : ℂ) := by
            rw [mul_comm ((starRingEnd ℂ) (PZ z w i)) (PZ z w i), Complex.mul_conj]
            push_cast; ring
  rw [stepA, stepB, Complex.ofReal_sum]
  exact Finset.sum_congr rfl (fun i _ => stepC i)

/-- **Bochner (atomic).** For nonnegative weights, the Toeplitz Hermitian form
    is real and `≥ 0` — the mechanism yielding `T^{g_0}, T^{f_ext} ⪰ 0`. -/
theorem bochner_toeplitz_form_nonneg
    (mu : Fin n → ℝ) (z : Fin n → ℂ) (hmu : ∀ i, 0 ≤ mu i) (w : Fin (K + 1) → ℂ) :
    0 ≤ (∑ j : Fin (K + 1), ∑ l : Fin (K + 1),
          (starRingEnd ℂ) (w j) * ToeplitzZ mu z j l * w l).re := by
  rw [bochner_toeplitz_form_eq, Complex.ofReal_re]
  exact Finset.sum_nonneg (fun i _ => mul_nonneg (hmu i) (Complex.normSq_nonneg _))

end
end ErdosV7
