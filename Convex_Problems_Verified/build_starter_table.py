"""Build a small starter cell table for the 6 multi-cell paper rows.
Each region is split into an N x N uniform grid in (h, p) with q held at
the row's full range. Output CSV is compatible with verify_all_emop.py.

Default grid is 4 x 4 (16 cells per region = 96 cells total).

Usage:
  python build_starter_table.py                       # 4x4 -> bnb/starter_4x4.csv
  python build_starter_table.py --grid 5              # 5x5 -> bnb/starter_5x5.csv
  python build_starter_table.py --grid 4 --out my.csv
"""
from __future__ import annotations
import argparse, csv, os, numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))

# Region definitions (h, p, q ranges + N, T, R config)
# A6/A7/A9 cover the same h/p box at different q strips; B0/B1/B2 are the
# narrow-q corner near p=0.4.
REGIONS = [
    {"row": "A6", "h_lo": 0.00, "h_hi": 0.08, "p_lo": 0.00, "p_hi": 1.00,
     "q_lo": -1.000, "q_hi": -0.050, "N": 10000, "T": 2000, "R": 10},
    {"row": "A7", "h_lo": 0.00, "h_hi": 0.08, "p_lo": 0.00, "p_hi": 1.00,
     "q_lo": -0.050, "q_hi": -0.025, "N": 10000, "T": 2000, "R": 10},
    {"row": "A9", "h_lo": 0.00, "h_hi": 0.08, "p_lo": 0.00, "p_hi": 1.00,
     "q_lo":  0.025, "q_hi":  0.050, "N": 10000, "T": 2000, "R": 10},
    {"row": "B0", "h_lo": 0.00, "h_hi": 0.06, "p_lo": 0.33, "p_hi": 0.45,
     "q_lo": -0.020, "q_hi":  0.020, "N": 10000, "T": 2000, "R": 10},
    {"row": "B1", "h_lo": 0.00, "h_hi": 0.06, "p_lo": 0.33, "p_hi": 0.45,
     "q_lo": -0.025, "q_hi": -0.020, "N": 10000, "T": 2000, "R": 10},
    {"row": "B2", "h_lo": 0.00, "h_hi": 0.06, "p_lo": 0.33, "p_hi": 0.45,
     "q_lo":  0.020, "q_hi":  0.025, "N": 10000, "T": 2000, "R": 10},
]


def grid_cells(region, n):
    """Tile region into n x n uniform cells in (h, p)."""
    h_edges = np.linspace(region["h_lo"], region["h_hi"], n + 1)
    p_edges = np.linspace(region["p_lo"], region["p_hi"], n + 1)
    out = []
    for i in range(n):
        for j in range(n):
            out.append({
                "row":   region["row"],
                "h1":    float(h_edges[i]),
                "h2":    float(h_edges[i + 1]),
                "p1":    float(p_edges[j]),
                "p2":    float(p_edges[j + 1]),
                "q1":    region["q_lo"],
                "q2":    region["q_hi"],
                "N":     region["N"],
                "T":     region["T"],
                "R":     region["R"],
                "level": "L0",
                "omega": "",       # to be filled by verifier
                "status":"pending",
                "src":   "starter",
            })
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--grid", type=int, default=4,
                    help="N for the N x N starter grid (default 4).")
    ap.add_argument("--out", default=None,
                    help="Output CSV path. Default bnb/starter_{N}x{N}.csv")
    args = ap.parse_args()

    n = args.grid
    out = args.out or os.path.join(HERE, "bnb", f"starter_{n}x{n}.csv")

    cells = []
    for region in REGIONS:
        cells.extend(grid_cells(region, n))

    header = ["row","h1","h2","p1","p2","q1","q2","N","T","R","level","omega","status","src"]
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=header)
        w.writeheader()
        for c in cells:
            w.writerow(c)

    print(f"Wrote {len(cells)} cells to {out}", flush=True)
    print(f"Breakdown:", flush=True)
    from collections import Counter
    for row, n_cells in Counter(c["row"] for c in cells).items():
        print(f"  {row}: {n_cells} cells", flush=True)


if __name__ == "__main__":
    main()
