#!/usr/bin/env python3
"""Table S1: baseline characteristics for the four TriNetX statin comparisons
(each statin vs a propensity-matched no-statin control; the analyses used in
the four-statin forest plot, Figure 2A).

Table S1A shows the cohorts before matching, Table S1B the matched cohorts.
Values are % of cohort (n in the column header) or mean ± SD; post-matching
SMDs >= 0.10 are bolded.

Step 04a of the pipeline.

Input:  data/raw/trinetx/<statin>_vs_no_statin/Baseline_Patient_Characteristics.csv
Output: results/tables/Table_S1_four_statins_vs_no_statin.{docx,csv}

Requires python-docx. Run from the repository root:
    python scripts/04a_table_s1_four_statins.py
"""

import csv
import io
from pathlib import Path

from docx import Document
from docx.enum.section import WD_ORIENT
from docx.enum.table import WD_CELL_VERTICAL_ALIGNMENT, WD_TABLE_ALIGNMENT
from docx.enum.text import WD_ALIGN_PARAGRAPH, WD_BREAK
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Inches, Pt, Twips

TRINETX_DIR = Path("data", "raw", "trinetx")
COMPARISONS = {
    statin: TRINETX_DIR / f"{statin.lower()}_vs_no_statin"
    for statin in ("Atorvastatin", "Rosuvastatin", "Simvastatin", "Pravastatin")
}
CONTROL_LABEL = "No statin"
OUT_DIR = Path("results", "tables")
OUT_DOCX = OUT_DIR / "Table_S1_four_statins_vs_no_statin.docx"
OUT_CSV = OUT_DIR / "Table_S1_four_statins_vs_no_statin.csv"

FONT = "Times New Roman"
BODY_PT = 8
SMD_FLAG = 0.10
LABEL_W, VALUE_W, SMD_W = 2900, 1000, 800  # twips; 2900 + 4 * 2800 = 14,100 (9.8 in)
COL_WIDTHS = [LABEL_W] + [VALUE_W, VALUE_W, SMD_W] * len(COMPARISONS)

# (label, characteristic ID, continuous?)
SECTIONS = [
    ("Demographics", [
        ("Age at index, years, mean ± SD", "AI", True),
        ("Male, %", "M", False),
        ("Female, %", "F", False),
    ]),
    ("Diagnoses, %", [
        ("Genitourinary diseases (N00–N99)", "N00-N99", False),
        ("Hypertensive diseases (I10–I1A)", "I10-I1A", False),
        ("Essential hypertension (I10)", "I10", False),
        ("Diabetes mellitus (E08–E13)", "E08-E13", False),
        ("Ischemic heart diseases (I20–I25)", "I20-I25", False),
        ("Other nervous system disorders (G89–G99)", "G89-G99", False),
        ("Heart failure (I50)", "I50", False),
        ("Neoplasms (C00–D49)", "C00-D49", False),
        ("Kidney transplant failure (T86.12)", "T86.12", False),
        ("Diseases of liver (K70–K77)", "K70-K77", False),
        ("Kidney transplant rejection (T86.11)", "T86.11", False),
        ("Tobacco use (Z72.0)", "Z72.0", False),
    ]),
    ("Medications, %", [
        ("Immunological agents (IM000)", "IM000", False),
        ("Beta-blocking agents (C07)", "C07", False),
        ("Aspirin", "1191", False),
        ("Renin–angiotensin agents (C09)", "C09", False),
    ]),
    ("Laboratory", [
        ("Total cholesterol recorded, %", "9000", False),
        ("Total cholesterol, mg/dL, mean ± SDᵃ", "9000", True),
    ]),
]


# ---------------------------------------------------------------------------
# Parse the TriNetX exports (preamble lines, then a plain CSV table)
# ---------------------------------------------------------------------------
def read_trinetx_table(path, header_prefix):
    lines = path.read_text(encoding="utf-8-sig").splitlines()
    start = next(i for i, line in enumerate(lines) if line.startswith(header_prefix))
    return list(csv.DictReader(io.StringIO("\n".join(lines[start:]))))


def find(records, char_id, continuous):
    """Continuous rows carry a mean; the duplicate row with the same ID gives
    the share of patients with a recorded value."""
    for r in records:
        if r["Characteristic ID"] == char_id and (r["Cohort 1 Before: Mean"] != "") == continuous:
            return r
    raise KeyError(f"{char_id} (continuous={continuous})")


data = {}
for statin, directory in COMPARISONS.items():
    records = read_trinetx_table(directory / "Baseline_Patient_Characteristics.csv", "Characteristic ID,")
    follow_up = read_trinetx_table(directory / "Follow_up_time_metrics_matched_results_table.csv", "Cohort,")
    age = find(records, "AI", True)
    data[statin] = {
        "records": records,
        "n": {(c, p): int(age[f"Cohort {c} {p}: Patient Count"]) for c in (1, 2) for p in ("Before", "After")},
        "median_follow_up": [r["Median Follow-up (Days)"] for r in follow_up],
    }


def fmt_value(r, cohort, phase, continuous):
    if continuous:
        mean = r[f"Cohort {cohort} {phase}: Mean"]
        if mean == "":
            return "—"
        return f"{float(mean):.1f} ± {float(r[f'Cohort {cohort} {phase}: SD']):.1f}"
    return f"{float(r[f'Cohort {cohort} {phase}: % of Cohort'].rstrip('%')):.1f}"


def smd_value(r, phase):
    x = r[f"{phase}: Standardized Mean Difference"]
    return None if x == "" else float(x)


# ---------------------------------------------------------------------------
# Word helpers
# ---------------------------------------------------------------------------
def set_cell_borders(cell, top=None, bottom=None):
    """Borders in eighths of a point; None means no border."""
    tc_pr = cell._tc.get_or_add_tcPr()
    borders = OxmlElement("w:tcBorders")
    for edge, size in (("top", top), ("left", None), ("bottom", bottom), ("right", None)):
        el = OxmlElement(f"w:{edge}")
        if size:
            el.set(qn("w:val"), "single")
            el.set(qn("w:sz"), str(size))
            el.set(qn("w:color"), "000000")
        else:
            el.set(qn("w:val"), "nil")
        borders.append(el)
    tc_pr.append(borders)


def shade_cell(cell, fill):
    shd = OxmlElement("w:shd")
    shd.set(qn("w:val"), "clear")
    shd.set(qn("w:color"), "auto")
    shd.set(qn("w:fill"), fill)
    cell._tc.get_or_add_tcPr().append(shd)


def write_cell(cell, text, bold=False, italic=False, left=False):
    cell.vertical_alignment = WD_CELL_VERTICAL_ALIGNMENT.CENTER
    para = cell.paragraphs[0]
    para.alignment = WD_ALIGN_PARAGRAPH.LEFT if left else WD_ALIGN_PARAGRAPH.CENTER
    para.paragraph_format.space_before = Pt(1)
    para.paragraph_format.space_after = Pt(1)
    run = para.add_run(text)
    run.bold = bold
    run.italic = italic
    run.font.size = Pt(BODY_PT)
    run.font.name = FONT


def add_paragraph(doc, runs, size, space_before=0, space_after=4):
    para = doc.add_paragraph()
    para.paragraph_format.space_before = Pt(space_before)
    para.paragraph_format.space_after = Pt(space_after)
    for text, bold in runs:
        run = para.add_run(text)
        run.bold = bold
        run.font.size = Pt(size)
        run.font.name = FONT
    return para


def fix_table_layout(table):
    tbl_pr = table._tbl.tblPr
    layout = OxmlElement("w:tblLayout")
    layout.set(qn("w:type"), "fixed")
    tbl_pr.append(layout)
    tbl_w = tbl_pr.find(qn("w:tblW"))
    tbl_w.set(qn("w:w"), str(sum(COL_WIDTHS)))
    tbl_w.set(qn("w:type"), "dxa")
    for grid_col, width in zip(table._tbl.tblGrid.findall(qn("w:gridCol")), COL_WIDTHS):
        grid_col.set(qn("w:w"), str(width))
    for row in table.rows:
        for i, width in enumerate(COL_WIDTHS):
            row.cells[i].width = Twips(width)


def mark_header_row(row):
    el = OxmlElement("w:tblHeader")
    el.set(qn("w:val"), "true")
    row._tr.get_or_add_trPr().append(el)


THICK, THIN = 12, 6
N_COLS = len(COL_WIDTHS)


def build_table(doc, phase, csv_rows):
    n_body = 1 + sum(len(rows) + 1 for _, rows in SECTIONS)  # +1 follow-up / blank row
    table = doc.add_table(rows=3 + n_body, cols=N_COLS)
    table.alignment = WD_TABLE_ALIGNMENT.CENTER
    table.autofit = False
    fix_table_layout(table)

    # Header row 1: statin comparison group labels
    hdr = table.rows[0]
    write_cell(hdr.cells[0], "")
    set_cell_borders(hdr.cells[0], top=THICK)
    for s, statin in enumerate(COMPARISONS):
        start = 1 + 3 * s
        merged = hdr.cells[start].merge(hdr.cells[start + 2])
        write_cell(merged, f"{statin} vs {CONTROL_LABEL.lower()}", bold=True)
        set_cell_borders(merged, top=THICK, bottom=THIN)

    # Header row 2: column names; row 3: cohort sizes
    write_cell(table.rows[1].cells[0], "Characteristic", bold=True, left=True)
    set_cell_borders(table.rows[1].cells[0])
    write_cell(table.rows[2].cells[0], "")
    set_cell_borders(table.rows[2].cells[0], bottom=THIN)
    for s, statin in enumerate(COMPARISONS):
        start = 1 + 3 * s
        n = data[statin]["n"]
        for k, (label, n_label) in enumerate([
            (statin, f"n = {n[(1, phase)]:,}"),
            (CONTROL_LABEL, f"n = {n[(2, phase)]:,}"),
            ("SMD", ""),
        ]):
            write_cell(table.rows[1].cells[start + k], label, bold=True)
            set_cell_borders(table.rows[1].cells[start + k])
            write_cell(table.rows[2].cells[start + k], n_label, italic=True)
            set_cell_borders(table.rows[2].cells[start + k], bottom=THIN)
    for row in table.rows[:3]:
        mark_header_row(row)

    r_idx = 3
    for title, rows in SECTIONS:
        merged = table.rows[r_idx].cells[0].merge(table.rows[r_idx].cells[N_COLS - 1])
        write_cell(merged, title, bold=True, left=True)
        set_cell_borders(merged)
        shade_cell(merged, "F2F2F2")
        r_idx += 1
        for label, char_id, continuous in rows:
            cells = table.rows[r_idx].cells
            write_cell(cells[0], "  " + label, left=True)
            set_cell_borders(cells[0])
            csv_row = [phase, title, label]
            for s, statin in enumerate(COMPARISONS):
                r = find(data[statin]["records"], char_id, continuous)
                smd = smd_value(r, phase)
                values = [fmt_value(r, 1, phase, continuous), fmt_value(r, 2, phase, continuous),
                          "—" if smd is None else f"{smd:.2f}"]
                flag = phase == "After" and smd is not None and round(smd, 2) >= SMD_FLAG
                for k, value in enumerate(values):
                    write_cell(cells[1 + 3 * s + k], value, bold=(k == 2 and flag))
                    set_cell_borders(cells[1 + 3 * s + k])
                csv_row += [values[0], values[1], "" if smd is None else f"{smd:.4f}"]
            csv_rows.append(csv_row)
            r_idx += 1

    # Final row: median follow-up (matched cohorts) or blank rule row
    cells = table.rows[r_idx].cells
    label = "Median follow-up, days" if phase == "After" else ""
    write_cell(cells[0], label, left=True)
    set_cell_borders(cells[0], top=THIN, bottom=THICK)
    for s, statin in enumerate(COMPARISONS):
        fu = data[statin]["median_follow_up"] if phase == "After" else ["", ""]
        for k, value in enumerate([fu[0], fu[1], ""]):
            write_cell(cells[1 + 3 * s + k], value)
            set_cell_borders(cells[1 + 3 * s + k], top=THIN, bottom=THICK)
    return table


# ---------------------------------------------------------------------------
# Build the document
# ---------------------------------------------------------------------------
doc = Document()
doc.styles["Normal"].font.name = FONT
doc.styles["Normal"].font.size = Pt(10)
section = doc.sections[0]
section.orientation = WD_ORIENT.LANDSCAPE
section.page_width, section.page_height = Inches(11), Inches(8.5)
section.left_margin = section.right_margin = Inches(0.6)
section.top_margin = section.bottom_margin = Inches(0.6)

footnote = (
    "Values are % of cohort unless otherwise indicated; cohort sizes are given in the column headers. "
    "Each statin cohort was compared with a separately propensity-matched (1:1) no-statin control "
    "cohort. SMD, absolute standardized mean difference; bold indicates post-matching SMD ≥ 0.10. "
    "Diagnoses are ICD-10-CM codes; medications are ATC classes except aspirin (RxNorm). "
    "ᵃAmong patients with a recorded value; mean cholesterol remained lower in statin users after "
    "matching. BMI percentile (recorded in ≤ 211 patients per cohort) and "
    "unknown sex (≤ 10 patients) are not shown."
)

csv_rows = []
for phase, suffix, caption in [
    ("Before", "A", "before propensity-score matching"),
    ("After", "B", "after 1:1 propensity-score matching"),
]:
    if phase == "After":
        doc.add_paragraph().add_run().add_break(WD_BREAK.PAGE)
    add_paragraph(doc, [
        (f"Table S1{suffix}. ", True),
        ("Baseline characteristics of kidney transplant recipients for each statin compared with "
         f"no-statin controls, {caption} (TriNetX US Collaborative Network).", False),
    ], size=11, space_after=6)
    build_table(doc, phase, csv_rows)
    add_paragraph(doc, [(footnote, False)], size=7.5, space_before=4, space_after=0)

OUT_DIR.mkdir(parents=True, exist_ok=True)
doc.save(OUT_DOCX)

with OUT_CSV.open("w", newline="", encoding="utf-8") as fh:
    writer = csv.writer(fh)
    header = ["Phase", "Section", "Characteristic"]
    for statin in COMPARISONS:
        header += [f"{statin}", f"{CONTROL_LABEL} ({statin} match)", f"SMD ({statin})"]
    writer.writerow(header)
    writer.writerows(csv_rows)

print(f"Wrote {OUT_DOCX} and {OUT_CSV}")
for statin, d in data.items():
    n = d["n"]
    after_smds = [smd_value(r, "After") for r in d["records"]
                  if r["Characteristic ID"] != "9084" and smd_value(r, "After") is not None]
    print(f"{statin:13s} before {n[(1, 'Before')]:>5,} vs {n[(2, 'Before')]:>6,}; "
          f"after {n[(1, 'After')]:>5,} per group; max post-matching SMD {max(after_smds):.3f}")
