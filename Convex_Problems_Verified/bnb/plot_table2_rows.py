"""Plot the per-row (h, p) verification grids for the 5 Table 2 rows.
One PNG per row, each saved as bnb/progress_<label>.png.

Cell coloring:
  - green: omega >= 0.38 (verified)
  - bright green: omega = inf (primal infeasible; vacuously holds)
  - red: omega < 0.38 (would-be failure; should be split or bumped)
  - gray hatched: still pending
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Patch

HERE = Path(__file__).resolve().parent
STATE_DIR = HERE / "table2_rows_state"

THRESHOLD = 0.3795

# Per-row (h, p, q) ranges for plot bounds and title
ROW_META = {
    "A6":   {"h": (0.00, 0.08), "p": (0.00, 1.00), "q": (-1.000, -0.050)},
    "A6_zoom": {"h": (0.008, 0.040), "p": (0.00, 0.40), "q": (-1.000, -0.050)},
    "A6_zoom_R50": {"h": (0.008, 0.040), "p": (0.00, 0.40), "q": (-1.000, -0.050)},
    "A6_zoom_R100": {"h": (0.008, 0.040), "p": (0.00, 0.40), "q": (-1.000, -0.050)},
    "A6_zoom_split_R50": {"h": (0.008, 0.040), "p": (0.00, 0.40), "q": (-1.000, -0.050)},
    "A6_1": {"h": (0.00, 0.08), "p": (0.00, 1.00), "q": (-1.000, -0.525)},
    "A6_2": {"h": (0.00, 0.08), "p": (0.00, 1.00), "q": (-0.525, -0.050)},
    "A7":   {"h": (0.00, 0.08), "p": (0.00, 1.00), "q": (-0.050, -0.025)},
    "A9":   {"h": (0.00, 0.08), "p": (0.00, 1.00), "q": ( 0.025,  0.050)},
    "B1":   {"h": (0.00, 0.06), "p": (0.33, 0.45), "q": (-0.025, -0.020)},
    "B2":   {"h": (0.00, 0.06), "p": (0.33, 0.45), "q": ( 0.020,  0.025)},
}


def plot_row(label: str):
    state_path = STATE_DIR / f"{label}.json"
    if not state_path.exists():
        print(f"  {label}: state file not yet present, skipping")
        return None

    boxes = json.load(open(state_path))["boxes"]
    leaves = [b for b in boxes if b["status"] != "split"]
    meta = ROW_META[label]

    fig, ax = plt.subplots(figsize=(12, 8))

    n_pass = n_inf = n_fail = n_pending = n_bumping = n_solver_err = 0
    for b in leaves:
        omega = b.get("omega")
        status = b["status"]
        is_inf = (omega is not None) and not math.isfinite(omega)
        if status == "pending" and omega is None:
            face, edge, lbl_color = "lightgray", "black", "black"
            hatch = "//"
            n_pending += 1
            text = f"L{b['level']}\n(pending)"
        elif status == "pending" and omega is not None:
            face, edge, lbl_color, hatch = "#ffd27f", "darkorange", "black", "\\\\"
            n_bumping += 1
            text = f"L{b['level']}\n{omega:.5f}"
        elif status == "verified" and is_inf:
            face, edge, lbl_color, hatch = "#9be39b", "darkgreen", "black", None
            n_inf += 1
            text = "inf"
        elif status == "verified":
            face, edge, lbl_color, hatch = "#7fcf7f", "darkgreen", "black", None
            n_pass += 1
            text = f"{omega:.5f}"
        elif status == "failed" and omega is None:
            # Solver crash / no result (NOT a threshold failure)
            face, edge, lbl_color, hatch = "#d8b3ff", "purple", "black", "xx"
            n_solver_err += 1
            text = "SOLVER\nFAIL"
        else:
            # status == "failed" with a finite omega < threshold
            face, edge, lbl_color, hatch = "#ff7f7f", "darkred", "black", None
            n_fail += 1
            text = f"{omega:.5f}"

        rect = Rectangle(
            (b["h_lo"], b["p_lo"]),
            b["h_hi"] - b["h_lo"], b["p_hi"] - b["p_lo"],
            facecolor=face, edgecolor=edge, linewidth=0.5, hatch=hatch,
        )
        ax.add_patch(rect)
        cx = 0.5 * (b["h_lo"] + b["h_hi"])
        cy = 0.5 * (b["p_lo"] + b["p_hi"])
        ax.text(cx, cy, text, ha="center", va="center",
                fontsize=5.5, color=lbl_color)

    ax.set_xlim(meta["h"])
    ax.set_ylim(meta["p"])
    ax.set_xlabel("h  (each cell is a band [h_lo, h_hi] enforced as h_1 ≤ E(M) ≤ h_2)")
    ax.set_ylabel("p  (each cell is a band [p_lo, p_hi] enforced as p_1 ≤ c_1 ≤ p_2)")
    total = len(leaves)
    n_verified = n_pass + n_inf
    title = (f"Row {label}: q = ({meta['q'][0]}, {meta['q'][1]}),  "
             f"threshold Ω ≥ {THRESHOLD}\n"
             f"verified {n_verified}/{total}  —  "
             f"Ω≥{THRESHOLD}: {n_pass}, infeasible: {n_inf}, "
             f"split/bump: {n_bumping}, fail: {n_fail}, "
             f"solver-err: {n_solver_err}, pending: {n_pending}")
    ax.set_title(title, fontsize=11)

    handles = [
        Patch(facecolor="#7fcf7f", edgecolor="darkgreen",
              label=f"Ω ≥ {THRESHOLD} (verified)"),
        Patch(facecolor="#9be39b", edgecolor="darkgreen",
              label="primal infeasible (vacuously verified)"),
        Patch(facecolor="#ffd27f", edgecolor="darkorange", hatch="\\\\",
              label="bump/split (Ω < threshold at lower N)"),
        Patch(facecolor="#ff7f7f", edgecolor="darkred",
              label="threshold failure (Ω < threshold)"),
        Patch(facecolor="#d8b3ff", edgecolor="purple", hatch="xx",
              label="solver failure (no result)"),
        Patch(facecolor="lightgray", edgecolor="black", hatch="//",
              label="pending (not solved yet)"),
    ]
    ax.legend(handles=handles, loc="lower right", framealpha=0.92, fontsize=8)

    plt.tight_layout()
    out = HERE / f"progress_{label}.png"
    plt.savefig(out, dpi=140)
    plt.close()
    print(f"  saved {out}  ({n_verified}/{total} verified, "
          f"{n_fail} fail, {n_pending} pending)")
    return out


def main():
    print(f"Plotting per-row progress to {HERE}/progress_<label>.png")
    for label in ["A6", "A6_zoom", "A6_zoom_R50", "A6_zoom_R100", "A6_zoom_split_R50", "A6_1", "A6_2", "A7", "A9", "B1", "B2"]:
        plot_row(label)


if __name__ == "__main__":
    main()
