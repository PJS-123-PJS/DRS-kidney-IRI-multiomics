#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Build Supplementary_Tables_S1-S6_20261006.docx (A4 landscape, header repeat,
monospace for primer sequence columns) from the six small-table CSVs."""
import csv
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
]

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
    for i, row in enumerate(data):
        for j, val in enumerate(row):
            cell = t.cell(i + 1, j)
            cell.text = ""
            r = cell.paragraphs[0].add_run(val)
            r.font.size = Pt(7.5)
            if j in mono_idx:
                r.font.name = "Consolas"
                r.element.rPr.rFonts.set(qn("w:eastAsia"), "Consolas")

out = f"{BASE}/Supplementary_Tables_S1-S6_20261006.docx"
doc.save(out)
print("saved:", out)
