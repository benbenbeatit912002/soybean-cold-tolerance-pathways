"""Read the author's local Table S2 without modifying it.

Requires openpyxl. The output is a comparison reference, never an analysis input.
Publisher identity has not been verified: the public download returned HTML.
"""
import argparse
import hashlib
import json
from pathlib import Path
import re
import zipfile
import xml.etree.ElementTree as ET
import openpyxl


def extract(path):
    workbook = openpyxl.load_workbook(path, read_only=True, data_only=True)
    records = []
    for sheet in workbook:
        method = sheet.title.split("(")[0].lower()
        phase = "short" if "Short" in sheet.title else "middle"
        for cells in sheet.iter_rows():
            # All eight source tables use B:G; do not infer rows from formatting.
            values = [c.value for c in cells]
            if len(values) < 7 or not re.fullmatch(r"GO:\d{7}", str(values[1])):
                continue
            if not all(isinstance(values[i], (int, float)) for i in (3, 4, 5, 6)):
                raise ValueError(f"Missing numeric reference at {sheet.title}!B{cells[1].row}:G{cells[1].row}")
            records.append(dict(phase=phase, method=method, go_id=values[1], description=values[2],
                annotation_size=values[3], ct_count=values[4], p_value=values[5], p_bh=values[6],
                source_sheet=sheet.title, source_range=f"B{cells[1].row}:G{cells[1].row}"))
    keys = [(r["phase"], r["method"], r["go_id"]) for r in records]
    if len(keys) != len(set(keys)):
        raise ValueError("Duplicate GO IDs in a source table")
    return {"source": str(path.resolve()), "source_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        "source_status": "Local manuscript supplement; publisher file identity unverified", "records": records}


def extract_selected(path, enrichment):
    """Read Table S1's printed probabilities, preserving less-than qualifiers."""
    ns = {"w": "http://schemas.openxmlformats.org/wordprocessingml/2006/main"}
    with zipfile.ZipFile(path) as archive:
        root = ET.fromstring(archive.read("word/document.xml"))
    normalize = lambda text: re.sub(r"[^a-z0-9]", "", text.lower())
    lookup = {(r["phase"], normalize(r["description"])): r["go_id"] for r in enrichment}
    def probability(text):
        value = text.replace(" ", "").replace("\u2212", "-")
        value = value.replace("\u00d710", "e").replace("<", "")
        return float(value)
    records = []
    for phase, table in zip(("short", "middle"), root.findall(".//w:tbl", ns)):
        for row_number, row in enumerate(table.findall("w:tr", ns)[2:], start=3):
            cells = ["".join(t.text or "" for t in c.findall(".//w:t", ns)) for c in row.findall("w:tc", ns)]
            if len(cells) < 11 or not cells[1].isdigit():
                continue
            go = lookup[(phase, normalize(cells[0]))]
            for method, pcol, qcol in (("hypergeometric",3,4), ("sumstat",6,7), ("gsea",9,10)):
                records.append(dict(phase=phase,method=method,go_id=go,description=cells[0],
                    annotation_size=int(cells[1]),ct_count=int(cells[2]),
                    p_value=probability(cells[pcol]),p_bh=probability(cells[qcol]),
                    p_less_than="<" in cells[pcol],q_less_than="<" in cells[qcol],
                    source_table=1 if phase=="short" else 2,source_row=row_number))
    return dict(source=str(path.resolve()),source_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),records=records)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("workbook", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--selected-docx", type=Path)
    args = parser.parse_args()
    data = extract(args.workbook)
    if args.selected_docx:
        data["selected"] = extract_selected(args.selected_docx, data["records"])
    args.output.write_text(json.dumps(data, indent=2, ensure_ascii=False)+"\n", encoding="utf-8")
    from collections import Counter
    print(dict(Counter((r["phase"], r["method"]) for r in data["records"])))
