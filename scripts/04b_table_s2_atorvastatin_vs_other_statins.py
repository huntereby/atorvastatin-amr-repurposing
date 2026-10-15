#!/usr/bin/env python3
"""Table S2: baseline characteristics of the atorvastatin and other-statin
TriNetX cohorts before and after 1:1 propensity-score matching.

Step 04b of the pipeline.

Input:  data/raw/trinetx/atorvastatin_vs_other_statins/Baseline_Patient_Characteristics.csv
Output: results/tables/Table_S2_atorvastatin_vs_other_statins.docx
        results/tables/Table_S2_atorvastatin_vs_other_statins.csv

Requires python-docx. Run from the repository root:
    python scripts/04b_table_s2_atorvastatin_vs_other_statins.py
"""

import csv
import io
import math
from pathlib import Path

from docx import Document
from docx.enum.section import WD_ORIENT
from docx.enum.table import WD_CELL_VERTICAL_ALIGNMENT, WD_TABLE_ALIGNMENT
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Inches, Pt, Twips

STUDY_DIR = Path("data", "raw", "trinetx", "atorvastatin_vs_other_statins")
INPUT_FILE = STUDY_DIR / "Baseline_Patient_Characteristics.csv"
FOLLOW_UP_FILE = STUDY_DIR / "Follow_up_time_metrics_matched_results_table.csv"
OUT_DIR = Path("results", "tables")
OUT_DOCX = OUT_DIR / "Table_S2_atorvastatin_vs_other_statins.docx"
OUT_CSV = OUT_DIR / "Table_S2_atorvastatin_vs_other_statins.csv"

FONT = "Times New Roman"
BODY_PT = 9
COL_WIDTHS = [4300, 1650, 1650, 900, 1650, 1650, 900]  # twips; 12,700 = 8.82 in


# ---------------------------------------------------------------------------
# Parse the TriNetX exports (preamble lines, then a plain CSV table)
# ---------------------------------------------------------------------------
def read_trinetx_table(path, header_prefix):
    lines = path.read_text(encoding="utf-8-sig").splitlines()
    start = next(i for i, line in enumerate(lines) if line.startswith(header_prefix))
    return list(csv.DictReader(io.StringIO("\n".join(lines[start:]))))


records = read_trinetx_table(INPUT_FILE, "Characteristic ID,")
follow_up = read_trinetx_table(FOLLOW_UP_FILE, "Cohort,")


def find(char_id, continuous=False):
    """Continuous rows carry a mean; the duplicate row with the same ID gives
    the number of patients with a recorded value."""
    for r in records:
        has_mean = r["Cohort 1 Before: Mean"] != ""
        if r["Characteristic ID"] == char_id and has_mean == continuous:
            return r
    raise KeyError(f"{char_id} (continuous={continuous})")


def fmt_int(x):
    return f"{int(float(x)):,}"


def fmt_pct(x):
    return f"{float(str(x).rstrip('%')):.1f}"


def fmt_smd(x):
    return "—" if x == "" else f"{float(x):.3f}"


def count(r, cohort, phase):
    return (f"{fmt_int(r[f'Cohort {cohort} {phase}: Patient Count'])} "
            f"({fmt_pct(r[f'Cohort {cohort} {phase}: % of Cohort'])})")


def mean_sd(r, cohort, phase):
    return (f"{float(r[f'Cohort {cohort} {phase}: Mean']):.1f} ± "
            f"{float(r[f'Cohort {cohort} {phase}: SD']):.1f}")


def row_values(r, fmt):
    return [
        fmt(r, 1, "Before"), fmt(r, 2, "Before"), fmt_smd(r["Before: Standardized Mean Difference"]),
        fmt(r, 1, "After"), fmt(r, 2, "After"), fmt_smd(r["After: Standardized Mean Difference"]),
    ]


age = find("AI", continuous=True)
cohort_n = [age[f"Cohort {c} {p}: Patient Count"] for p in ("Before", "After") for c in (1, 2)]

SECTIONS = [
    ("Demographics", [
        ("Age at index, years, mean ± SD", age, mean_sd),
        ("Male, n (%)", find("M"), count),
        ("Female, n (%)", find("F"), count),
    ]),
    ("Diagnoses, n (%)", [
        ("Diseases of the genitourinary system (N00–N99)", find("N00-N99"), count),
        ("Hypertensive diseases (I10–I1A)", find("I10-I1A"), count),
        ("Essential (primary) hypertension (I10)", find("I10"), count),
        ("Diabetes mellitus (E08–E13)", find("E08-E13"), count),
        ("Ischemic heart diseases (I20–I25)", find("I20-I25"), count),
        ("Other disorders of the nervous system (G89–G99)", find("G89-G99"), count),
        ("Heart failure (I50)", find("I50"), count),
        ("Neoplasms (C00–D49)", find("C00-D49"), count),
        ("Kidney transplant failure (T86.12)", find("T86.12"), count),
        ("Diseases of liver (K70–K77)", find("K70-K77"), count),
        ("Kidney transplant rejection (T86.11)", find("T86.11"), count),
        ("Tobacco use (Z72.0)", find("Z72.0"), count),
    ]),
    ("Medications, n (%)", [
        ("Immunological agents (IM000)", find("IM000"), count),
        ("Beta-blocking agents (C07)", find("C07"), count),
        ("Aspirin", find("1191"), count),
        ("Agents acting on the renin–angiotensin system (C09)", find("C09"), count),
    ]),
    ("Laboratory", [
        ("Total cholesterol recorded, n (%)", find("9000"), count),
        ("Total cholesterol, mg/dL, mean ± SDᵃ", find("9000", continuous=True), mean_sd),
    ]),
]

max_after_smd = max(
    float(r["After: Standardized Mean Difference"])
    for r in records if r["After: Standardized Mean Difference"] != ""
)
median_follow_up = sorted({r["Median Follow-up (Days)"] for r in follow_up})


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


def mark_header_row(row):
    tr_pr = row._tr.get_or_add_trPr()
    el = OxmlElement("w:tblHeader")
    el.set(qn("w:val"), "true")
    tr_pr.append(el)


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
section.top_margin = section.bottom_margin = Inches(0.75)

add_paragraph(doc, [
    ("Table S2. ", True),
    ("Baseline characteristics of kidney transplant recipients receiving atorvastatin or other "
     "statins, before and after 1:1 propensity-score matching (TriNetX US Collaborative Network).", False),
], size=11, space_after=8)

n_body = sum(len(rows) + 1 for _, rows in SECTIONS)
table = doc.add_table(rows=3 + n_body, cols=7)
table.alignment = WD_TABLE_ALIGNMENT.CENTER
table.autofit = False

# Fixed layout with explicit grid and table widths so Word keeps the columns
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

THICK, THIN = 12, 6

# Header row 1: spanning group labels
hdr = table.rows[0]
write_cell(hdr.cells[0], "")
set_cell_borders(hdr.cells[0], top=THICK)
for start, label in ((1, "Before matching"), (4, "After matching")):
    merged = hdr.cells[start].merge(hdr.cells[start + 2])
    write_cell(merged, label, bold=True)
    set_cell_borders(merged, top=THICK, bottom=THIN)

# Header row 2: column names
for i, label in enumerate(["Characteristic", "Atorvastatin", "Other statins", "SMD",
                           "Atorvastatin", "Other statins", "SMD"]):
    write_cell(table.rows[1].cells[i], label, bold=True, left=(i == 0))
    set_cell_borders(table.rows[1].cells[i])

# Header row 3: cohort sizes
n_labels = ["", f"n = {fmt_int(cohort_n[0])}", f"n = {fmt_int(cohort_n[1])}", "",
            f"n = {fmt_int(cohort_n[2])}", f"n = {fmt_int(cohort_n[3])}", ""]
for i, label in enumerate(n_labels):
    write_cell(table.rows[2].cells[i], label, italic=True)
    set_cell_borders(table.rows[2].cells[i], bottom=THIN)

for row in table.rows[:3]:
    mark_header_row(row)

csv_rows = []
r_idx = 3
for s_idx, (title, rows) in enumerate(SECTIONS):
    merged = table.rows[r_idx].cells[0].merge(table.rows[r_idx].cells[6])
    write_cell(merged, title, bold=True, left=True)
    set_cell_borders(merged)
    shade_cell(merged, "F2F2F2")
    r_idx += 1
    for j, (label, record, fmt) in enumerate(rows):
        last = s_idx == len(SECTIONS) - 1 and j == len(rows) - 1
        values = row_values(record, fmt)
        cells = table.rows[r_idx].cells
        write_cell(cells[0], "    " + label, left=True)
        for k, value in enumerate(values, start=1):
            write_cell(cells[k], value)
        for c in cells:
            set_cell_borders(c, bottom=THICK if last else None)
        csv_rows.append([title, label, *values])
        r_idx += 1

smd_bound = math.ceil(max_after_smd * 100) / 100
cholesterol_smd = float(find("9000", continuous=True)["After: Standardized Mean Difference"])
cholesterol_note = (
    f"; mean cholesterol remained lower in atorvastatin users after matching (SMD {cholesterol_smd:.2f})"
    if cholesterol_smd >= 0.10 else ""
)
add_paragraph(doc, [(
    "Values are n (% of cohort) unless otherwise indicated. SMD, absolute standardized mean "
    f"difference; all post-matching SMDs were < {smd_bound:.2f}. Diagnoses are ICD-10-CM codes; "
    "medications are ATC classes except aspirin (RxNorm).", False)], size=8, space_before=6, space_after=2)
add_paragraph(doc, [(
    "ᵃAmong patients with a recorded value" + cholesterol_note + ". BMI percentile (recorded in 10 atorvastatin patients and "
    "no comparators before matching, and in no patients after matching) and unknown sex (0 patients) "
    f"are not shown. Median follow-up was {' / '.join(median_follow_up)} days in both matched cohorts.",
    False)], size=8, space_after=0)

OUT_DIR.mkdir(parents=True, exist_ok=True)
doc.save(OUT_DOCX)

with OUT_CSV.open("w", newline="", encoding="utf-8") as fh:
    writer = csv.writer(fh)
    writer.writerow(["Section", "Characteristic",
                     f"Atorvastatin before (n={cohort_n[0]})", f"Other statins before (n={cohort_n[1]})",
                     "SMD before",
                     f"Atorvastatin after (n={cohort_n[2]})", f"Other statins after (n={cohort_n[3]})",
                     "SMD after"])
    writer.writerows(csv_rows)

print(f"Wrote {OUT_DOCX} and {OUT_CSV}")
print(f"Cohorts before: {cohort_n[0]} vs {cohort_n[1]}; after: {cohort_n[2]} vs {cohort_n[3]}")
print(f"Max post-matching SMD: {max_after_smd:.4f}; median follow-up: {median_follow_up}")
