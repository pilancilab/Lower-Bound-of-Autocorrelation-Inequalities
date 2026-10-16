"""Plot the 100 A6_zoom cells colored by which (N, T, R) config first
verified them, with Ω printed in each cell.

The 4-tier lineage:
  driver L0   (N=10000,  T=2000,  R=10)   — cells verified at original solve
  R=50 retry  (N=10000,  T=2000,  R=50)   — cells failing L0 but cleared at R=50
  R=100 retry (N=10000,  T=2000,  R=100)  — (empty: helped no cell)
  NT + R=30   (N=30000,  T=6000,  R=30)   — cells needing the larger discretization

Plus 1 "still under 0.3795" cell (cleared the relaxed 0.3794 bound).
"""
from __future__ import annotations
import json
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Patch

ROOT = Path(__file__).resolve().parent
STATE = ROOT / "table2_rows_state"

THRESHOLD_STRICT  = 0.3795
THRESHOLD_RELAXED = 0.3794

# (N, T, R) configs per provenance stage
CONFIGS = {
    "driver_L0":   (10000,  2000, 10),
    "R50":         (10000,  2000, 50),
    "R100":        (10000,  2000, 100),
    "NT_R30":      (30000,  6000, 30),
    "still_below": (30000,  6000, 30),   # same config, still < 0.3795
}

# Color per stage (light → dark = more compute needed)
COLORS = {
    "driver_L0":   "#a8e6a3",   # light green
    "R50":         "#ffe28a",   # yellow
    "R100":        "#ffb060",   # orange (unused)
    "NT_R30":      "#7fbcff",   # blue
    "still_below": "#ff7f7f",   # red
}

STAGE_FILES = [
    ("driver_L0", "A6_zoom.json"),
    ("R50",       "A6_zoom_R50.json"),
    ("R100",      "A6_zoom_R100.json"),
    ("NT_R30",    "A6_zoom_NT_R30.json"),
]


def main():
    # Load state files
    stages = {name: json.load(open(STATE/fname))["boxes"]
              for name, fname in STAGE_FILES}

    # Final boxes (from A6_zoom_NT_R30) — the truth for Ω
    final = {(b["h_lo"], b["h_hi"], b["p_lo"], b["p_hi"]): b
             for b in stages["NT_R30"]}

    # Provenance: which stage first marked each cell as "verified"
    provenance = {}
    for key, b_final in final.items():
        for name, _ in STAGE_FILES:
            b = next((x for x in stages[name]
                      if (x["h_lo"], x["h_hi"], x["p_lo"], x["p_hi"]) == key), None)
            if b and b.get("status") == "verified":
                provenance[key] = name
                break
        else:
            # Never verified at any stage; check if it's the relaxed-pass cell
            if (b_final.get("omega") is not None
                    and b_final["omega"] >= THRESHOLD_RELAXED):
                provenance[key] = "still_below"
            else:
                provenance[key] = "still_below"  # treat as red regardless

    # Plot
    fig, ax = plt.subplots(figsize=(13, 10))
    for key, b in final.items():
        stage = provenance[key]
        face = COLORS[stage]
        rect = Rectangle(
            (b["h_lo"], b["p_lo"]),
            b["h_hi"] - b["h_lo"],
            b["p_hi"] - b["p_lo"],
            facecolor=face, edgecolor="black", linewidth=0.6,
        )
        ax.add_patch(rect)
        cx = 0.5 * (b["h_lo"] + b["h_hi"])
        cy = 0.5 * (b["p_lo"] + b["p_hi"])
        omega = b.get("omega")
        omega_str = f"{omega:.5f}" if omega is not None else "—"
        ax.text(cx, cy, omega_str, ha="center", va="center",
                fontsize=7.5, color="black")

    ax.set_xlim(0.008, 0.04)
    ax.set_ylim(0.0, 0.4)
    ax.set_xlabel(r"$h$  (cell-uniform $h_1 = h_2$)")
    ax.set_ylabel(r"$p$  (cell-uniform $p_1 = p_2$)")

    n_by_stage = {s: 0 for s in COLORS}
    for s in provenance.values():
        n_by_stage[s] += 1

    title = ("A6_zoom certification (q ∈ [-1.000, -0.050])\n"
             f"100 cells; global min Ω = "
             f"{min(b['omega'] for b in final.values() if b.get('omega') is not None):.6f}; "
             f"certified Ω ≥ 0.3794")
    ax.set_title(title)

    def lbl(name, count):
        N, T, R = CONFIGS[name]
        if name == "still_below":
            return f"Ω<0.3795 (cleared 0.3794): {count}"
        return f"{name} (N={N}, T={T}, R={R}): {count}"

    handles = [Patch(facecolor=COLORS[s], edgecolor="black", label=lbl(s, n_by_stage[s]))
               for s in ["driver_L0", "R50", "R100", "NT_R30", "still_below"]
               if n_by_stage[s] > 0 or s == "R100"]
    ax.legend(handles=handles, loc="upper right", framealpha=0.95, fontsize=9)

    plt.tight_layout()
    out = ROOT / "progress_A6_zoom_NT_R30.png"
    plt.savefig(out, dpi=140)
    print(f"saved {out}")
    print(f"cells by stage: {n_by_stage}")


if __name__ == "__main__":
    main()
