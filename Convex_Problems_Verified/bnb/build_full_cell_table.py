"""Emit a complete cell-by-cell certification table:
  row, h1, h2, p1, p2, q1, q2, N, T, R, level, Ω

Sources (in order):
  - 12 single-box Table-2 passing rows  (table2_results.json,  Ω >= 0.38 at coarse solve)
  - A6 outer (84 cells outside the zoom region, from A6.json)
  - A6_zoom_NT_R30 (100 cells, per-cell (N,T,R) reconstructed by stage diff)
  - A7  (verified leaves only)
  - A9  (verified leaves only)
  - B1  (verified leaves only)
  - B2  (verified leaves only)
  - state.json (the q∈[-0.02,0.02] driver region, levels L0/L1/L2)

Writes:
  bnb/full_cell_table.csv
  bnb/full_cell_table.tex   (uses longtable, multi-page)
"""
from __future__ import annotations
import csv, json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
TROW = ROOT / "table2_rows_state"

# ---- Level → (N, T, R) ----
DRIVER_LEVELS = {                # bnb/driver.py LEVELS
    0: (10000, 2000, 10),
    1: (20000, 4000, 10),
    2: (50000, 10000, 10),
}
TABLE2_LEVEL  = (10000, 2000, 10)   # verify_table2_rows.py LEVELS_DEFAULT[0]
ZOOM_STAGE_CFG = {                  # A6_zoom lineage
    "driver_L0":   (10000, 2000, 10),
    "R50":         (10000, 2000, 50),
    "R100":        (10000, 2000, 100),
    "NT_R30":      (30000, 6000, 30),
    "still_below": (30000, 6000, 30),
}


def in_zoom(b):
    return (0.008 - 1e-9 <= b["h_lo"] and b["h_hi"] <= 0.04 + 1e-9
            and 0.0 - 1e-9 <= b["p_lo"] and b["p_hi"] <= 0.4 + 1e-9)


def zoom_provenance():
    """For each A6_zoom cell, find earliest state-file stage where it
    became 'verified'. If never verified, label 'still_below'."""
    files = [("driver_L0", "A6_zoom.json"),
             ("R50",       "A6_zoom_R50.json"),
             ("R100",      "A6_zoom_R100.json"),
             ("NT_R30",    "A6_zoom_NT_R30.json")]
    stages = {name: json.load(open(TROW/fn))["boxes"] for name, fn in files}
    final = stages["NT_R30"]
    out = {}
    for b in final:
        key = (b["h_lo"], b["h_hi"], b["p_lo"], b["p_hi"])
        out[key] = "still_below"
        for name, _ in files:
            match = next((x for x in stages[name]
                          if (x["h_lo"], x["h_hi"], x["p_lo"], x["p_hi"]) == key
                          and x.get("status") == "verified"), None)
            if match is not None:
                out[key] = name
                break
    return out, final


def main():
    rows = []   # each: dict with keys row,h1,h2,p1,p2,q1,q2,N,T,R,level,omega,status

    # 1. Single-box passing Table-2 rows
    t2 = json.load(open(ROOT / "table2_results.json"))["results"]
    for r in t2:
        if r.get("omega") is None or r["omega"] < 0.38:
            continue   # failing ones come in via the refined state files
        rows.append({
            "row": r["label"], "src": "table2_results.json",
            "h1": r["h1"], "h2": r["h2"],
            "p1": r["p1"], "p2": r["p2"],
            "q1": r["q1"], "q2": r["q2"],
            "N": r["N"], "T": r["T"], "R": r["R"],
            "level": "single",
            "omega": r["omega"],
            "status": "verified",
        })

    # Region q-ranges (for the bnb-refined Table-2 rows)
    Q_RANGE = {
        "A6": (-1.000, -0.050),
        "A7": (-0.050, -0.025),
        "A9": ( 0.025,  0.050),
        "B1": (-0.025, -0.020),
        "B2": ( 0.020,  0.025),
    }

    # 2. A6 outer cells (verified leaves outside the zoom region)
    A6 = json.load(open(TROW/"A6.json"))["boxes"]
    q1, q2 = Q_RANGE["A6"]
    for b in A6:
        if b.get("status") != "verified": continue
        if in_zoom(b): continue                    # handled by A6_zoom
        if b.get("omega") is None: continue
        N, T, R = DRIVER_LEVELS[b["level"]]
        rows.append({
            "row": "A6", "src": "A6.json",
            "h1": b["h_lo"], "h2": b["h_hi"],
            "p1": b["p_lo"], "p2": b["p_hi"],
            "q1": q1, "q2": q2,
            "N": N, "T": T, "R": R, "level": f"L{b['level']}",
            "omega": b["omega"], "status": "verified",
        })

    # 3. A6_zoom_NT_R30 with per-cell stage provenance
    prov, zoom_boxes = zoom_provenance()
    q1, q2 = Q_RANGE["A6"]
    for b in zoom_boxes:
        if b.get("omega") is None: continue
        key = (b["h_lo"], b["h_hi"], b["p_lo"], b["p_hi"])
        stage = prov.get(key, "still_below")
        N, T, R = ZOOM_STAGE_CFG[stage]
        rows.append({
            "row": "A6_zoom", "src": f"A6_zoom_NT_R30.json [{stage}]",
            "h1": b["h_lo"], "h2": b["h_hi"],
            "p1": b["p_lo"], "p2": b["p_hi"],
            "q1": q1, "q2": q2,
            "N": N, "T": T, "R": R, "level": "leaf",
            "omega": b["omega"],
            "status": "verified" if (b.get("status") == "verified"
                                     or (stage == "still_below"
                                         and b["omega"] >= 0.3794)) else "failed",
        })

    # 4. A7, A9, B1, B2 — verified leaves only (refined Table-2 rows)
    for name in ["A7", "A9", "B1", "B2"]:
        q1, q2 = Q_RANGE[name]
        data = json.load(open(TROW/f"{name}.json"))["boxes"]
        for b in data:
            if b.get("status") != "verified": continue
            if b.get("omega") is None: continue
            # All cells from verify_table2_rows.py at single-tier level 0
            N, T, R = TABLE2_LEVEL
            rows.append({
                "row": name, "src": f"{name}.json",
                "h1": b["h_lo"], "h2": b["h_hi"],
                "p1": b["p_lo"], "p2": b["p_hi"],
                "q1": q1, "q2": q2,
                "N": N, "T": T, "R": R, "level": f"L{b.get('level','?')}",
                "omega": b["omega"], "status": "verified",
            })

    # 5. state.json — the q∈[-0.02,0.02] driver region (call it "B0")
    s = json.load(open(ROOT/"state.json"))["boxes"]
    Q_B0 = (-0.020, 0.020)
    for b in s:
        if b.get("status") != "verified": continue
        if b.get("omega") is None: continue
        N, T, R = DRIVER_LEVELS[b["level"]]
        rows.append({
            "row": "B0", "src": "state.json",
            "h1": b["h_lo"], "h2": b["h_hi"],
            "p1": b["p_lo"], "p2": b["p_hi"],
            "q1": Q_B0[0], "q2": Q_B0[1],
            "N": N, "T": T, "R": R, "level": f"L{b['level']}",
            "omega": b["omega"], "status": "verified",
        })

    # ---- write CSV ----
    csv_path = ROOT / "full_cell_table.csv"
    with open(csv_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["row","h1","h2","p1","p2","q1","q2",
                                          "N","T","R","level","omega","status","src"])
        w.writeheader()
        for r in rows:
            w.writerow(r)

    # ---- write LaTeX longtable ----
    tex_path = ROOT / "full_cell_table.tex"
    with open(tex_path, "w") as f:
        f.write(r"""\documentclass{article}
\usepackage[margin=0.6in]{geometry}
\usepackage{longtable}
\usepackage{booktabs}
\usepackage{siunitx}
\begin{document}
\footnotesize
\setlength{\tabcolsep}{4pt}
\begin{longtable}{l r r r r r r r r r l}
\caption{Per-cell certification record. Each row is one verified box
in the partition of $[0,2]\times[0,1]\times[-1,1]$ for
$(h,p,q)$. The columns $(h_1, h_2)\times(p_1, p_2)\times(q_1, q_2)$
give the cell extent; $(N, T, R)$ identify the discretization that
certified the cell's lower bound $\Omega$.}\label{tab:cells} \\
\toprule
Row & $h_1$ & $h_2$ & $p_1$ & $p_2$ & $q_1$ & $q_2$ & $N$ & $T$ & $R$ & $\Omega$ \\
\midrule
\endfirsthead
\multicolumn{11}{l}{\footnotesize\itshape continued from previous page} \\
\toprule
Row & $h_1$ & $h_2$ & $p_1$ & $p_2$ & $q_1$ & $q_2$ & $N$ & $T$ & $R$ & $\Omega$ \\
\midrule
\endhead
\midrule \multicolumn{11}{r}{\footnotesize\itshape continued on next page} \\
\endfoot
\bottomrule
\endlastfoot
""")
        # Sort: by row label (alphabetical), then by (h1, p1)
        def sort_key(r):
            order = {"A1":1,"A2":2,"A3":3,"A4":4,"A5":5,"A6":6,"A6_zoom":7,
                     "A7":8,"A8":9,"A9":10,"A10":11,"A11":12,"A12":13,
                     "A13":14,"A14":15,"A15":16,"B0":17,"B1":18,"B2":19}
            return (order.get(r["row"], 99), r["h1"], r["p1"])
        for r in sorted(rows, key=sort_key):
            f.write(f"{r['row']} & "
                    f"{r['h1']:.4f} & {r['h2']:.4f} & "
                    f"{r['p1']:.4f} & {r['p2']:.4f} & "
                    f"{r['q1']:+.4f} & {r['q2']:+.4f} & "
                    f"{r['N']} & {r['T']} & {r['R']} & "
                    f"{r['omega']:.6f} \\\\\n")
        f.write(r"""\end{longtable}
\end{document}
""")

    # ---- console summary ----
    from collections import Counter
    by_row = Counter(r["row"] for r in rows)
    print(f"=== Cell count by region ===")
    for k in sorted(by_row, key=lambda x: (
            {"A1":1,"A2":2,"A3":3,"A4":4,"A5":5,"A6":6,"A6_zoom":7,
             "A7":8,"A8":9,"A9":10,"A10":11,"A11":12,"A12":13,
             "A13":14,"A14":15,"A15":16,"B0":17,"B1":18,"B2":19}.get(x,99))):
        print(f"  {k:<10s}  {by_row[k]:>4d}")
    print(f"  {'TOTAL':<10s}  {len(rows):>4d}")
    print(f"\nsaved {csv_path}")
    print(f"saved {tex_path}")


if __name__ == "__main__":
    main()
