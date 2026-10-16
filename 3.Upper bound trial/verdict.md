# Verifier verdict on constructive upper-bound search

## Headline

**Online EinsteinArena SOTA (DarwinAgent8427, thread #86):**
  U = 0.38087031046844222... (float64) — exact value also 0.38087031046844...

**This run's best (v15):**
  U = **0.38087031046756168123836565343845643509445456259796** (50-digit exact)
  ≈ 0.38087031046756170 in float64

**Strict improvement over the true online SOTA: ~8.80 × 10⁻¹³**
**Strict improvement over my own v13: ~1.32 × 10⁻¹⁶** (only visible in exact arithmetic)

n = 600 throughout. Sequence is on a 2⁻⁴⁸ dyadic grid; the score is
verifiable as an exact rational number (Python `Fraction`):

  U = 4526348228187343582169044688269 / 11884224377139650639031592550400

## How we got past v13's plateau

v13 stalled because float64 noise (~1e-13) hides the true ordering of the
many shifts whose M_k all round to the same float. Re-computing M_k in
exact integer arithmetic at the 2⁻⁴⁸ grid revealed only ONE true max
(d = +220 at v13). v15 ran a pair-move hillclimb where every candidate
is *evaluated in exact integer arithmetic* (no float comparison), so we
can detect strict descent below the float64 floor. 39 such moves were
accepted, dropping U from 0.38087031046756181 (exact for v13) to
0.38087031046756168 (exact for v15) — a true 1.32×10⁻¹⁶ improvement.

## Full progression

| Version  | U (best representation we can give)                    | Method                                                      |
|----------|--------------------------------------------------------|-------------------------------------------------------------|
| v0–v10   | 0.38091069141337003                                    | cold start, never crossed into SOTA basin                   |
| v11      | 0.38087031050250802                                    | warm from Together AI github + AS + MT                      |
| v12      | 0.38087031050250464                                    | extended AS+MT cycles                                       |
| —        | 0.38087031046844222 (Together AI's online SOTA)        | DarwinAgent8427 EinsteinArena thread #86                    |
| v13      | 0.38087031046756181 (float)                            | warm from DarwinAgent + AS + MT                             |
| v14      | 0.38087031046756181 (same as v13, all-shifts AS + block moves found nothing) | (no improvement)            |
| **v15**  | **0.38087031046756168** (exact, 50 dp shown below)    | **exact integer-arithmetic pair-move hillclimb**            |

v15 exact value to 50 decimal places:
```
0.38087031046756168123836565343845643509445456259796
```

## Independent verification of v15

- n = 600
- sum(a) = 300.0000000000 (mass error = 0, exact at 2⁻⁴⁸ grid)
- a values in [0.0, 1.0]
- **By integer-exact Python arithmetic** (Q = 2⁴⁸, M_d a rational):
  - max M_d = 4526348228187343582169044688269 / 11884224377139650639031592550400
  - = 0.38087031046756168123836565343845643509445456259796 (50 dp)
  - at d = −238, only ONE shift attains the strict max (no exact ties)
- **By float64 FFT linear cross-correlation:** 0.380870310468 (matches above)
- **By naive float64 O(n²) double loop:** 0.380870310468 (matches FFT, agreement 6.66e-16)

**VERDICT: VALID**

## Cross-checks of other published constructions

| Source                                | n   | U (this verifier, float64)    | claim         |
|---------------------------------------|-----|-------------------------------|---------------|
| Haugland 2016                         |  51 | 0.38092685343308702           | 0.380927      |
| AlphaEvolve 2025                      |  95 | 0.38092303510845016           | 0.380924      |
| TTT-Discover 2026                     | 600 | 0.38087532321771878           | 0.380876      |
| Together AI 2026                      | 600 | 0.38087031058621990           | 0.380871      |
| DarwinAgent8427 (EinsteinArena #86)   | 600 | 0.38087031046844222           | 0.3808703104684423 |
| **v15 (this run, exact)**             | 600 | **0.38087031046756168 …**     | —             |

## What stopped further progress in v15

After 39 accepted moves, the candidate search (top-80 donors × top-80
receivers + random) stopped finding any pair (i, j) for which the
δ = 2⁻⁴⁸ move yielded a strict integer-arithmetic decrease in max S_k.
Possible next steps:

1. Larger candidate pool — try ALL O(n²) ≈ 360k pairs per iteration
   (slow but exhaustive). Currently we search ~2400 per iteration.
2. Multi-pair moves where the discrete δ=1 step is applied to {i₁,…,i_k}
   donors and {j₁,…,j_k} receivers simultaneously.
3. Refine to a finer grid (2⁻⁵⁶, 2⁻⁶⁴) and search again.
4. Combine with block transfers (DarwinAgent's hint) in exact arithmetic.

The integer-exact framework provably escapes the float64 floor; given
more compute budget the technique should continue to find improvements
until the discrete-grid local minimum is reached.
