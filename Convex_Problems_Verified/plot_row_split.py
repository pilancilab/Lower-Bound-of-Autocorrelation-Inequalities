"""Plot a row with 4x4 parents + nested 2x2 splits for failing cells.
Goes as deep as we have data for (regions_4x4, regions_split_2x2, regions_split_2x2x2)."""
import os, json, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Patch

HERE = "/home/sykim/Theory_Harness/Convex_Problems_Verified"
THRESHOLD = 0.37912

# Try to load each level (some files may not exist yet).
PARENT_JSONS = [
    os.path.join(HERE, "regions_4x4_results.json"),
    os.path.join(HERE, "regions_B0_4x4_results.json"),
    os.path.join(HERE, "regions_A6_4x4_results.json"),
]
LEVEL_JSONS = [
    os.path.join(HERE, "regions_B0_split_2x2_results.json"),
    os.path.join(HERE, "regions_B0_split_2x2x2_results.json"),
    os.path.join(HERE, "regions_B0_split_2x2x2x2_results.json"),
    os.path.join(HERE, "regions_A6_split_2x2_results.json"),
    os.path.join(HERE, "regions_A6_split_2x2x2_results.json"),
    os.path.join(HERE, "regions_A6_bumped_N30k_T6k_R30_results.json"),
    os.path.join(HERE, "regions_A6_split_after_bump_results.json"),
    os.path.join(HERE, "regions_split_2x2_results.json"),
    os.path.join(HERE, "regions_split_2x2x2_results.json"),
    os.path.join(HERE, "regions_split_2x2x2x2_results.json"),
    os.path.join(HERE, "regions_split_2x2x2x2x2_results.json"),
]

ROW_CFG = {
    "A6": (-1.000, -0.050, 0.08, 0.0, 1.00),
    "A7": (-0.050, -0.025, 0.08, 0.0, 1.00),
    "A9": ( 0.025,  0.050, 0.08, 0.0, 1.00),
    "B0": (-0.020,  0.020, 0.06, 0.33, 0.45),
    "B1": (-0.025, -0.020, 0.06, 0.33, 0.45),
    "B2": ( 0.020,  0.025, 0.06, 0.33, 0.45),
}

COLOR_VERIFIED   = "#5fa85f"
COLOR_INFEASIBLE = "#bfe3bf"
COLOR_BELOW_THR  = "#f08080"
COLOR_SOLVER_ERR = "#d8b4fe"
COLOR_PENDING    = "lightgray"


def load_all_levels(row):
    """Build a flat dict label -> record covering all levels for this row."""
    out = {}

    # Level 0: 4x4 parents (merged across all parent JSONs)
    for pjson in PARENT_JSONS:
        if not os.path.exists(pjson):
            continue
        for p in json.load(open(pjson)):
            lab = p.get("label", "")
            if lab.startswith(row + "_"):
                out[lab] = p

    # Levels 1, 2: split JSONs (and bumped-NTR JSONs, which use bump_label
    # like "A6_1_0.10.00#N30000T6000R30" — strip the "#N..." suffix so they
    # override the corresponding sub_label record from the prior split).
    for path in LEVEL_JSONS:
        if not os.path.exists(path):
            continue
        for s in json.load(open(path)):
            lab = s.get("sub_label") or s.get("bump_label", "")
            if "#" in lab:
                lab = lab.split("#", 1)[0]
            if lab.startswith(row + "_"):
                out[lab] = s
    return out


def status_of(record):
    """Return (status, value) where status in {CERT, BELOW, INF, ERR, PENDING}."""
    if record is None:
        return ("PENDING", None)
    L = record.get("certified_lower_bound")
    err = record.get("error") or ""
    if L is not None:
        if L == float("inf"):
            return ("INF", None)
        if L >= THRESHOLD:
            return ("CERT", L)
        return ("BELOW", L)
    if "zero-size array" in err:
        return ("INF", None)
    if err:
        return ("ERR", err)
    return ("PENDING", None)


COLORS = {
    "CERT":    (COLOR_VERIFIED,   None,  "black"),
    "INF":     (COLOR_INFEASIBLE, "inf", "black"),
    "BELOW":   (COLOR_BELOW_THR,  None,  "black"),
    "ERR":     (COLOR_SOLVER_ERR, "err", "black"),
    "PENDING": (COLOR_PENDING,    "pending", "dimgray"),
}


def draw_recursive(ax, label, h_l, h_h, p_l, p_h, records, depth=0):
    """Draw the cell `label`. If it's failing AND has 4 children in `records`,
    recurse on the children. Else draw a single block."""
    status, val = status_of(records.get(label))

    # If this cell is failing (BELOW or ERR), look for split children
    child_labels = [f"{label}.{ii}{jj}" for ii in range(2) for jj in range(2)]
    has_children = any(records.get(cl) is not None for cl in child_labels)

    if status in ("BELOW", "ERR") and has_children:
        # Recurse: thicker outline to mark the parent being split
        ax.add_patch(Rectangle((h_l, p_l), h_h - h_l, p_h - p_l,
                                facecolor="none", edgecolor="black", lw=1.5 + 0.5 * (3 - depth)))
        hmid = 0.5 * (h_l + h_h)
        pmid = 0.5 * (p_l + p_h)
        for ii, (sh_l, sh_h) in enumerate([(h_l, hmid), (hmid, h_h)]):
            for jj, (sp_l, sp_h) in enumerate([(p_l, pmid), (pmid, p_h)]):
                child = f"{label}.{ii}{jj}"
                draw_recursive(ax, child, sh_l, sh_h, sp_l, sp_h, records, depth + 1)
        return

    # Leaf: draw single block
    color, default_text, tc = COLORS[status]
    if default_text is None:
        text = f"{val:.5f}" if val is not None else "?"
    else:
        text = default_text
    # Smaller font for deeper levels
    fontsize = max(6, 10 - 2 * depth)
    ax.add_patch(Rectangle((h_l, p_l), h_h - h_l, p_h - p_l,
                            facecolor=color, edgecolor="white", lw=1.2))
    ax.text(0.5 * (h_l + h_h), 0.5 * (p_l + p_h), text,
            ha="center", va="center", color=tc, fontsize=fontsize)


def count_leaves(records, top_labels):
    """Count leaf statuses across the entire tree for the title."""
    counts = {"CERT": 0, "INF": 0, "BELOW": 0, "ERR": 0, "PENDING": 0}
    def visit(label):
        status, _ = status_of(records.get(label))
        child_labels = [f"{label}.{ii}{jj}" for ii in range(2) for jj in range(2)]
        has_children = any(records.get(cl) is not None for cl in child_labels)
        if status in ("BELOW", "ERR") and has_children:
            for c in child_labels:
                visit(c)
        else:
            counts[status] += 1
    for top in top_labels:
        visit(top)
    return counts


def plot(row, out_path):
    q_lo, q_hi, h_hi, p_lo, p_hi = ROW_CFG[row]
    h_edges = np.linspace(0.0, h_hi, 5)
    p_edges = np.linspace(p_lo, p_hi, 5)
    records = load_all_levels(row)

    fig, ax = plt.subplots(figsize=(12, 8))
    top_labels = []
    for i in range(4):
        for j in range(4):
            top_labels.append(f"{row}_{i}_{j}")
            draw_recursive(ax, f"{row}_{i}_{j}",
                           h_edges[i], h_edges[i + 1],
                           p_edges[j], p_edges[j + 1],
                           records, depth=0)

    counts = count_leaves(records, top_labels)

    ax.set_xlim(0.0, h_hi)
    ax.set_ylim(p_lo, p_hi)
    ax.set_xlabel("h  (each cell band [h_lo, h_hi])")
    ax.set_ylabel("p  (each cell band [p_lo, p_hi])")
    ax.set_title(f"Row {row}: q = ({q_lo:+.3f}, {q_hi:+.3f}),  threshold Ω ≥ {THRESHOLD}\n"
                 f"4×4 + adaptive 2×2 splits — "
                 f"cert={counts['CERT']}, inf={counts['INF']}, "
                 f"below={counts['BELOW']}, err={counts['ERR']}, "
                 f"pending={counts['PENDING']}")

    handles = [
        Patch(facecolor=COLOR_VERIFIED,   label=f"Ω ≥ {THRESHOLD} (verified)"),
        Patch(facecolor=COLOR_INFEASIBLE, label="primal infeasible (vacuously verified)"),
        Patch(facecolor=COLOR_BELOW_THR,  label=f"below threshold (Ω < {THRESHOLD})"),
        Patch(facecolor=COLOR_SOLVER_ERR, label="solver failure"),
        Patch(facecolor=COLOR_PENDING,    label="pending"),
    ]
    ax.legend(handles=handles, loc="lower right", framealpha=0.9, fontsize=9)
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    print(f"wrote {out_path}", flush=True)
    plt.close(fig)


if __name__ == "__main__":
    rows = sys.argv[1:] if len(sys.argv) > 1 else ["A7", "A9", "B1", "B2"]
    for r in rows:
        out = os.path.join(HERE, f"progress_{r}_split.png")
        plot(r, out)
