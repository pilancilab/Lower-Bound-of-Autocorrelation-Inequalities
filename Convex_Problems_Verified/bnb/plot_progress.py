"""Plot BnB progress: each box colored binary by Omega vs threshold,
with the Omega value printed inside."""
from __future__ import annotations
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Patch

HERE = Path(__file__).resolve().parent
STATE = HERE / "state.json"
OUT = HERE / "progress.png"

THRESHOLD = 0.3795


def main():
    boxes = json.load(open(STATE))["boxes"]
    leaves = [b for b in boxes if b["status"] != "split"]

    fig, ax = plt.subplots(figsize=(11, 8))

    for b in leaves:
        omega = b.get("omega")
        passed = (omega is not None) and (omega >= THRESHOLD)
        if omega is None:
            face, edge, text_color = "lightgray", "black", "black"
        elif passed:
            face, edge, text_color = "#7fcf7f", "darkgreen", "black"
        else:
            face, edge, text_color = "#ff7f7f", "darkred", "black"

        rect = Rectangle(
            (b["h_lo"], b["p_lo"]),
            b["h_hi"] - b["h_lo"],
            b["p_hi"] - b["p_lo"],
            facecolor=face, edgecolor=edge, linewidth=0.6,
        )
        ax.add_patch(rect)

        cx = 0.5 * (b["h_lo"] + b["h_hi"])
        cy = 0.5 * (b["p_lo"] + b["p_hi"])
        if omega is None:
            label = f"L{b['level']}\n(pending)"
        else:
            label = f"{omega:.5f}"
        ax.text(cx, cy, label, ha="center", va="center",
                fontsize=6.0, color=text_color)

    ax.set_xlim(0.0, 0.06)
    ax.set_ylim(0.33, 0.45)
    ax.set_xlabel("h  (h_1 = h_2 within each box)")
    ax.set_ylabel("p  (p_1 = p_2 within each box)")

    n_pass = sum(1 for b in leaves
                 if b.get("omega") is not None and b["omega"] >= THRESHOLD)
    n_fail = sum(1 for b in leaves
                 if b.get("omega") is not None and b["omega"] < THRESHOLD)
    n_pend = sum(1 for b in leaves if b.get("omega") is None)
    title = (f"EMOP BnB: Omega vs threshold={THRESHOLD}   |   "
             f"pass={n_pass}  fail={n_fail}  pending={n_pend}   |   "
             f"q=(-0.02, 0.02)")
    ax.set_title(title)

    handles = [
        Patch(facecolor="#7fcf7f", edgecolor="darkgreen",
              label=f"Omega >= {THRESHOLD}"),
        Patch(facecolor="#ff7f7f", edgecolor="darkred",
              label=f"Omega <  {THRESHOLD}"),
        Patch(facecolor="lightgray", edgecolor="black",
              label="pending (not yet solved)"),
    ]
    ax.legend(handles=handles, loc="lower right", framealpha=0.95, fontsize=8)

    plt.tight_layout()
    plt.savefig(OUT, dpi=140)
    print(f"saved {OUT}  ({n_pass} pass / {n_fail} fail / {n_pend} pending)")


if __name__ == "__main__":
    main()
