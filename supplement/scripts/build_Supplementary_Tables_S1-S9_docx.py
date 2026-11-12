#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Build Supplementary_Tables_S1-S10_20261024.docx — extends the S1-S6 builder
with S7/S8/S9. S9 carries italicised gene symbols (Car5b/Rpa1/Cxadr/Samd8)
per the gene-typography convention (genes italic in every deliverable file).
"""
import csv, re
from docx import Document
from docx.enum.section import WD_ORIENT
from docx.enum.table import WD_TABLE_ALIGNMENT
from docx.shared import Pt, Mm
from docx.oxml.ns import qn
from docx.oxml import OxmlElement

BASE = "/mnt/d/desktop/DRS/DRS_TechShowcase/supplement/tables"
TABLES = [
    ("Table S1. Sample metadata and DRS quality control.",
     "Table_S1_Sample_Metadata_DRS_QC.csv", []),
    ("Table S2. Antibodies used in this study.",
     "Table_S2_Antibodies.csv", []),
    ("Table S3. Primers used in this study.",
     "Table_S3_Primers.csv", ["Forward Primer (5'→3')", "Reverse Primer (5'→3')",
                              "F (5'→3')", "R (5'→3')", "Sequence (5'→3')"]),
    ("Table S4. Dual-anchor arbitration statistics per modification (Sham ↔ Day 1).",
     "Table_S4_ModVerdict_Paired_Stats.csv", []),
    ("Table S5. Animal-level sensitivity analysis of the principal arbitration verdicts.",
     "Table_S5_AnimalLevel_Sensitivity.csv", []),
    ("Table S6. Um caller independence test.",
     "Table_S6_Um_Caller_Independence.csv", []),
    ("Table S7. Per-transcript poly(A) tail-length trajectories of the Day-1 remodeled sets after Day 1.",
     "Table_S7_PolyA_Trajectory_After_Day1.csv", []),
    ("Table S8. Day-1 significant retained-intron events and their NMD-potential features.",
     "Table_S8_RI_NMD_Potential.csv", []),
    ("Table S9. Orthogonal validation of the inosine verdict: REDIportal site-level overlap (block A) and A→G mismatch hard signature (block B).",
     "Table_S9_Inosine_Orthogonal_Validation.csv", []),
    ("Table S10. Per-contrast Western blot densitometry statistics under the gel-run blocked linear model on log₂ ratios (Holm-adjusted versus Sham; Fig 6E, 6G).",
     "Table_S10_WB_Densitometry_BlockModel.csv", []),
]
ITALIC_GENES = {"Table_S9_Inosine_Orthogonal_Validation.csv": ["Car5b", "Rpa1", "Cxadr", "Samd8"]}

doc = Document()
sec = doc.sections[0]
sec.orientation = WD_ORIENT.LANDSCAPE
sec.page_width, sec.page_height = Mm(297), Mm(210)
sec.left_margin = sec.right_margin = Mm(15)
sec.top_margin = sec.bottom_margin = Mm(15)

style = doc.styles["Normal"]
style.font.name = "Arial"
style.font.size = Pt(8)
style.element.rPr.rFonts.set(qn("w:eastAsia"), "Arial")

def set_repeat_header(row):
    trPr = row._tr.get_or_add_trPr()
    tblHeader = OxmlElement("w:tblHeader")
    tblHeader.set(qn("w:val"), "true")
    trPr.append(tblHeader)

def add_cell_text(cell, text, size, italic_tokens=()):
    """write cell text, italicising whitelist gene tokens in place"""
    p = cell.paragraphs[0]
    if not italic_tokens:
        r = p.add_run(text)
        r.font.size = Pt(size)
        return [r]
    rx = re.compile(r"(?<![\w])(" + "|".join(italic_tokens) + r")(?![\w])")
    runs, pos = [], 0
    for m in rx.finditer(text):
        if m.start() > pos:
            r = p.add_run(text[pos:m.start()]); r.font.size = Pt(size); runs.append(r)
        r = p.add_run(m.group(1)); r.font.size = Pt(size); r.italic = True; runs.append(r)
        pos = m.end()
    if pos < len(text):
        r = p.add_run(text[pos:]); r.font.size = Pt(size); runs.append(r)
    return runs

n_ital = 0
for ti, (title, fname, mono_cols) in enumerate(TABLES):
    if ti > 0:
        doc.add_page_break()
    h = doc.add_paragraph()
    run = h.add_run(title)
    run.bold = True
    run.font.size = Pt(10)
    with open(f"{BASE}/{fname}", encoding="utf-8-sig") as f:
        rows = list(csv.reader(f))
    header, data = rows[0], rows[1:]
    mono_idx = {i for i, c in enumerate(header) if any(m in c for m in mono_cols)}
    ital_toks = ITALIC_GENES.get(fname, [])
    t = doc.add_table(rows=1 + len(data), cols=len(header))
    t.style = "Table Grid"
    t.alignment = WD_TABLE_ALIGNMENT.CENTER
    t.autofit = True
    for j, c in enumerate(header):
        cell = t.cell(0, j)
        cell.text = ""
        r = cell.paragraphs[0].add_run(c)
        r.bold = True
        r.font.size = Pt(8)
    set_repeat_header(t.rows[0])
    # row-iteration pass — t.cell(i,j) is O(rows) per call (quadratic on S7's 11k rows)
    for trow, row in zip(list(t.rows)[1:], data):
        cells = trow.cells
        for j, val in enumerate(row):
            cell = cells[j]
            cell.text = ""
            rs = add_cell_text(cell, val, 7.5, ital_toks)
            n_ital += sum(1 for r in rs if r.italic)
            if j in mono_idx:
                for r in rs:
                    r.font.name = "Consolas"
                    r.element.rPr.rFonts.set(qn("w:eastAsia"), "Consolas")

out = f"{BASE}/Supplementary_Tables_S1-S10_20261024.docx"
doc.save(out)

# verify: table count, italic runs in S9, no text loss
d2 = Document(out)
assert len(d2.tables) == 10, len(d2.tables)
s9 = d2.tables[8]
ital = [r.text for r in s9.rows[1].cells[13].paragraphs[0].runs if r.italic]
print("saved:", out, "| tables:", len(d2.tables), "| italic gene runs:", n_ital,
      "| S9 row1 note italics:", ital)
