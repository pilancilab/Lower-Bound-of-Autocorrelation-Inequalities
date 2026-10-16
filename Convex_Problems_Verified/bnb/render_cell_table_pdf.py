"""Render bnb/full_cell_table.csv to a multi-page PDF using fpdf2.
Fallback when no LaTeX compiler is available."""
from __future__ import annotations
import csv
from pathlib import Path
from fpdf import FPDF

ROOT = Path(__file__).resolve().parent
CSV  = ROOT / "full_cell_table.csv"
OUT  = ROOT / "full_cell_table.pdf"

# Sort order: same as the LaTeX builder
ORDER = {"A1":1,"A2":2,"A3":3,"A4":4,"A5":5,"A6":6,"A6_zoom":7,
         "A7":8,"A8":9,"A9":10,"A10":11,"A11":12,"A12":13,
         "A13":14,"A14":15,"A15":16,"B0":17,"B1":18,"B2":19}

# Column layout: (header, width_mm, format_fn)
COLS = [
    ("Row",   16, lambda r: r["row"]),
    ("h1",    16, lambda r: f"{float(r['h1']):.4f}"),
    ("h2",    16, lambda r: f"{float(r['h2']):.4f}"),
    ("p1",    16, lambda r: f"{float(r['p1']):.4f}"),
    ("p2",    16, lambda r: f"{float(r['p2']):.4f}"),
    ("q1",    18, lambda r: f"{float(r['q1']):+.4f}"),
    ("q2",    18, lambda r: f"{float(r['q2']):+.4f}"),
    ("N",     16, lambda r: r["N"]),
    ("T",     14, lambda r: r["T"]),
    ("R",     10, lambda r: r["R"]),
    ("Omega", 22, lambda r: f"{float(r['omega']):.6f}"),
]
ROW_H = 4.2     # mm
HDR_H = 5.5     # mm


class PDF(FPDF):
    def header(self):
        if self.page_no() == 1:
            self.set_font("Helvetica", "B", 11)
            self.cell(0, 6, "Per-cell certification record (748 cells)", new_x="LMARGIN", new_y="NEXT")
            self.set_font("Helvetica", "", 8)
            self.multi_cell(0, 4,
                "Each row is one verified box in the partition of (h,p,q). "
                "Columns (h1,h2)x(p1,p2)x(q1,q2) give the cell extent; "
                "(N,T,R) the discretization that certified the cell's lower bound Omega.")
            self.ln(1)
        self.set_font("Helvetica", "B", 7.5)
        self.set_fill_color(220, 220, 220)
        for hdr, w, _ in COLS:
            self.cell(w, HDR_H, hdr, border=1, align="C", fill=True)
        self.ln(HDR_H)
        self.set_font("Helvetica", "", 7)

    def footer(self):
        self.set_y(-10)
        self.set_font("Helvetica", "I", 7)
        self.cell(0, 4, f"Page {self.page_no()}", align="C")


def main():
    rows = []
    with open(CSV) as f:
        for r in csv.DictReader(f):
            rows.append(r)
    rows.sort(key=lambda r: (ORDER.get(r["row"], 99), float(r["h1"]), float(r["p1"])))

    pdf = PDF(orientation="L", unit="mm", format="A4")  # landscape for the wide table
    pdf.set_auto_page_break(auto=True, margin=12)
    pdf.add_page()

    for idx, r in enumerate(rows):
        # Color the one failing cell (Omega < 0.3795) in light red
        omega = float(r["omega"])
        if omega < 0.3795:
            pdf.set_fill_color(255, 215, 215); fill=True
        else:
            fill=False
        for hdr, w, fmt in COLS:
            pdf.cell(w, ROW_H, fmt(r), border=1, align="R", fill=fill)
        pdf.ln(ROW_H)

    pdf.output(str(OUT))
    print(f"saved {OUT}  ({len(rows)} rows, {pdf.page_no()} pages)")


if __name__ == "__main__":
    main()
