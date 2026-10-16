"""Faithful Python/CVXPY port of White's OPL program
(White-Minimum-Overlap/Programs/MO_convex_primal.txt).

Constraint-by-constraint translation, preserving:
  - paper-typo eps_bound (leading `2*m` with m = loop index)
  - L*L*sum(...)*N/2 form in (5.4)  [equivalent to L*sum]
  - NO j^2 lower-bound on second moment
  - NO Toeplitz / Bochner-Herglotz block
  - qsq_max = q1*q1 (uses q1^2, not max(q1^2, q2^2))

Inputs: same as EMOP.py.
"""
from __future__ import annotations

import os
os.environ.setdefault("OMP_NUM_THREADS", "8")
os.environ.setdefault("MKL_NUM_THREADS", "8")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "8")

import json, sys, time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT))

import numpy as np
import cvxpy as cp


def build_and_solve(N, T, R, h1, h2, p1, p2, q1, q2, threads=8):
    """Faithful port of MO_convex_primal.txt. Returns dict with omega and stats."""
    pi = np.pi
    L = 2.0 / N

    # ---- coefficient arrays (alphaminus, alphaplus, betaminus, betaplus) ----
    # 1-indexed via padding: arr[j][m] for j=1..N, m=1..2R
    j_arr = np.arange(1, N + 1)[:, None]       # (N, 1)
    m_arr = np.arange(1, 2 * R + 1)[None, :]   # (1, 2R)
    center = pi * m_arr * (j_arr - 0.5) / N
    half = pi * m_arr * 0.5 / N
    alphaminus = np.cos(center) - half  # (N, 2R)
    alphaplus  = np.cos(center) + half
    betaplus   = np.sin(center) + half
    betaminus  = np.sin(center) - half

    # ---- signvec[j] = (-1)^(j+1) for j=1..2T (1-indexed)  ----
    # In OPL, signvec[k+1] for k=1..T  =>  (-1)^k for k=1..T.
    # We just compute (-1)^k for k=1..T directly.
    sign_k = np.array([(-1.0) ** k for k in range(1, T + 1)])  # length T, index k=1..T
    # signvec[m] for m=1..R  =>  (-1)^(m+1)
    sign_m = np.array([(-1.0) ** (m + 1) for m in range(1, R + 1)])  # length R, index m=1..R

    # ---- epsilonbounds, deltabounds (paper-typo form) ----
    m_R = np.arange(1, R + 1, dtype=float)
    denom = 4.0 - (m_R / T) ** 2
    epsilonbounds = (1.0 / pi) * (1.0 / denom) * 2.0 * m_R * (6.0 * T ** 3) ** -0.5
    deltabounds   = (4.0 / pi) * (1.0 / denom) * (2.0 * T) ** -0.5

    # ---- decision variables ----
    Omega = cp.Variable(nonneg=True, name="Omega")
    w = cp.Variable(N, nonneg=True, name="w")  # w[j-1] = OPL's w[j], j=1..N
    v = cp.Variable(N, nonneg=True, name="v")
    c = cp.Variable(T, name="c")               # c[k-1] = OPL's c[k], k=1..T
    d = cp.Variable(T, name="d")
    epsilon = cp.Variable(R, name="epsilon")   # epsilon[m-1] = OPL's epsilon[m]
    delta   = cp.Variable(R, name="delta")

    # ---- derived expressions a_odd, a_even, b_odd, b_even (m=1..R, 0-indexed) ----
    constraints = []
    a_odd_list = []
    b_odd_list = []
    a_even_list = []
    b_even_list = []
    for m in range(1, R + 1):
        mm = 2 * m - 1
        # Σ_{k=1..T} signvec[k+1] * c[k] / (mm^2 - 4 k^2)
        denom_k = mm * mm - 4 * (np.arange(1, T + 1) ** 2)  # length T, may include zeros for mm=2k
        # Note: for mm odd and k integer, mm^2 - 4k^2 is never zero (mm odd, 4k^2 even diff). Safe.
        coef_c = sign_k / denom_k  # length T
        sum_c = coef_c @ c  # CVXPY expression
        a_odd_m = (epsilon[m - 1]
                   + 2 * mm * sign_m[m - 1] * (1.0 / pi)
                     * (0.5 / (mm * mm) + sum_c))
        a_odd_list.append(a_odd_m)

        coef_d = (np.arange(1, T + 1, dtype=float) * sign_k) / denom_k  # length T
        sum_d = coef_d @ d
        b_odd_m = (delta[m - 1] + (4.0 / pi) * sign_m[m - 1] * sum_d)
        b_odd_list.append(b_odd_m)

        a_even_list.append(0.5 * c[m - 1])   # OPL's a_even[m] = 0.5*c[m]
        b_even_list.append(0.5 * d[m - 1])

    # ---- box bounds on w_j, v_j ----
    constraints += [w <= Omega, v <= Omega]

    # ---- (5.2) L * Σ (w+v) == 1 ----
    constraints += [L * cp.sum(w + v) == 1]

    # ---- (5.3) L^2 * Σ (j*w - (j-1)*v) >= h1 ----
    j_vec = np.arange(1, N + 1, dtype=float)
    constraints += [L * L * (j_vec @ w - (j_vec - 1.0) @ v) >= h1]

    # ---- (5.4) variant: L * L * Σ (j-1)^2 (w+v) <= (2/3 + h2^2/2) * N/2 ----
    # Note: (N/2) on RHS is OPL's form; equivalent to L*Σ <= 2/3 + h2^2/2 since L=2/N.
    constraints += [
        L * L * ((j_vec - 1.0) ** 2 @ (w + v)) <= (2.0 / 3.0 + h2 * h2 / 2.0) * (N / 2.0)
    ]

    # ---- forall m in 1..R block ----
    for m in range(1, R + 1):
        idx_2m   = 2 * m - 1        # column index for alphaminus[:, 2m-1]: 0-indexed = 2m-1
        idx_2mm1 = 2 * m - 2        # column index for alphaminus[:, 2m-2]: 0-indexed = 2m-2  (this is OPL's 2m-1)
        # Wait: OPL alphaminus[j][2*m] -> 0-indexed column 2m-1
        #       OPL alphaminus[j][2*m-1] -> 0-indexed column 2m-2
        col_2m       = 2 * m - 1
        col_2m_minus = 2 * m - 2

        a_e = a_even_list[m - 1]
        a_o = a_odd_list[m - 1]
        b_e = b_even_list[m - 1]
        b_o = b_odd_list[m - 1]

        # 0.5*L * sum alphaminus[j][2m] * (w+v) <= -2*(a_e^2 + b_e^2)
        constraints += [
            0.5 * L * (alphaminus[:, col_2m] @ (w + v))
            <= -2.0 * (cp.square(a_e) + cp.square(b_e))
        ]
        # 0.5*L * sum alphaminus[j][2m-1] * (w+v) <= 4*signvec[m]*a_o/(pi*(2m-1)) - 2*(a_o^2 + b_o^2)
        constraints += [
            0.5 * L * (alphaminus[:, col_2m_minus] @ (w + v))
            <= 4.0 * sign_m[m - 1] * a_o / (pi * (2 * m - 1))
               - 2.0 * (cp.square(a_o) + cp.square(b_o))
        ]
        # betaminus[j][2m]*w - betaplus[j][2m]*v summed <= 0
        constraints += [
            0.5 * L * (betaminus[:, col_2m] @ w - betaplus[:, col_2m] @ v) <= 0
        ]
        # for (2m-1): RHS = -signvec[m]*b_o*8/(pi*(2m-1))
        constraints += [
            0.5 * L * (betaminus[:, col_2m_minus] @ w - betaplus[:, col_2m_minus] @ v)
            <= -sign_m[m - 1] * b_o * 8.0 / (pi * (2 * m - 1))
        ]
        # >= form
        constraints += [
            0.5 * L * (betaplus[:, col_2m] @ w - betaminus[:, col_2m] @ v) >= 0
        ]
        constraints += [
            0.5 * L * (betaplus[:, col_2m_minus] @ w - betaminus[:, col_2m_minus] @ v)
            >= -sign_m[m - 1] * b_o * 8.0 / (pi * (2 * m - 1))
        ]
        # epsilon, delta bounds
        constraints += [
            epsilon[m - 1] <=  epsilonbounds[m - 1],
            epsilon[m - 1] >= -epsilonbounds[m - 1],
            delta[m - 1]   <=  deltabounds[m - 1],
            delta[m - 1]   >= -deltabounds[m - 1],
        ]

    # ---- c, d bounds ----
    constraints += [c <= 2.0 / pi, c >= -2.0 / pi, d <= 2.0 / pi, d >= -2.0 / pi]
    constraints += [cp.sum_squares(c) + cp.sum_squares(d) <= 0.5]

    # ---- (5.12) p1 <= c1, c1 <= p2, q1 <= d1, d1 <= q2 ----
    constraints += [p1 <= c[0], c[0] <= p2, q1 <= d[0], d[0] <= q2]

    # ---- (5.13) qsq_max: OPL line 130 uses q1*q1 (NOT max(q1^2, q2^2)) ----
    constraints += [
        0.5 * L * (alphaplus[:, 1] @ (w + v))  # column index 1 = OPL's alphaplus[j][2]
        >= -0.5 * (p2 * p2 + q1 * q1)
    ]

    # ---- solve ----
    prob = cp.Problem(cp.Minimize(Omega), constraints)
    t0 = time.time()
    prob.solve(solver=cp.MOSEK, verbose=False,
               mosek_params={"MSK_IPAR_NUM_THREADS": int(threads)})
    elapsed = time.time() - t0

    # max violation
    max_violation = 0.0
    for con in constraints:
        try:
            vv = con.violation()
            if vv is None:
                continue
            mv = float(np.max(np.abs(vv)))
            max_violation = max(max_violation, mv)
        except Exception:
            pass

    return {
        "status": prob.status,
        "omega": prob.value,
        "max_violation": max_violation,
        "elapsed_s": elapsed,
    }


def solve_row(args):
    label, N, T, R, h1, h2, p1, p2, q1, q2, threads = args
    try:
        r = build_and_solve(N, T, R, h1, h2, p1, p2, q1, q2, threads)
        return {"label": label, "N": N, "T": T, "R": R, "ok": True, **r}
    except Exception as e:
        import traceback
        return {"label": label, "N": N, "T": T, "R": R, "ok": False,
                "error": f"{type(e).__name__}: {e}",
                "traceback": traceback.format_exc()}


def main():
    THRESHOLD = 0.38
    # 5 previously-failing rows at paper's resolution
    ROWS = [
        ("A6", 10000, 4000, 10, 0.00, 0.08, 0.00, 1.00, -1.000, -0.050),
        ("A7", 10000, 4000, 10, 0.00, 0.08, 0.00, 1.00, -0.050, -0.025),
        ("A9", 10000, 4000, 10, 0.00, 0.08, 0.00, 1.00,  0.025,  0.050),
        ("B1", 20000, 5000, 10, 0.00, 0.06, 0.33, 0.45, -0.025, -0.020),
        ("B2", 20000, 5000, 10, 0.00, 0.06, 0.33, 0.45,  0.020,  0.025),
    ]
    args = [(*r, 8) for r in ROWS]
    print(f"=== {len(args)} rows via faithful OPL port (4 workers x 8 threads) ===",
          flush=True)
    results = []
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=4) as ex:
        futs = {ex.submit(solve_row, a): a[0] for a in args}
        for fut in as_completed(futs):
            r = fut.result()
            results.append(r)
            if r.get("ok"):
                v = "PASS" if r["omega"] >= THRESHOLD else "FAIL"
                print(f"  {r['label']:>4s}  N={r['N']:>5d} T={r['T']:>5d}  "
                      f"omega={r['omega']:.6f}  mv={r['max_violation']:.1e}  "
                      f"{v}  ({r['elapsed_s']:.1f}s)", flush=True)
            else:
                print(f"  {r['label']:>4s}  ERROR {r.get('error')}", flush=True)
                print(r.get("traceback", "")[:600])
    print(f"Wall: {time.time()-t0:.1f}s")
    out = HERE / "white_opl_port_results.json"
    with open(out, "w") as f:
        json.dump({"threshold": THRESHOLD, "results": results}, f, indent=1)
    print(f"saved {out}")


if __name__ == "__main__":
    main()
