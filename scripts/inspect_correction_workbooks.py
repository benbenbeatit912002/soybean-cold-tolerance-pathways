"""Read-only formula archaeology. Does not modify or recalculate Excel files.

Requires openpyxl for this optional audit only; the score calculators do not.
"""
import argparse
import hashlib
import json
from pathlib import Path
import openpyxl

def inspect(path):
    workbook = openpyxl.load_workbook(path, data_only=False, read_only=True)
    result = {"path": str(path), "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
              "sheets": []}
    for sheet in workbook:
        record = {"sheet": sheet.title, "rows": sheet.max_row, "columns": sheet.max_column,
                  "formula_count": 0, "formula_examples": [], "target_rows": [], "header": []}
        for index, row in enumerate(sheet.iter_rows(), 1):
            if index <= 3:
                record["header"].append([str(c.value) for c in row[:15]])
            target = any("14g129500" in str(c.value).lower() or "14g041700" in str(c.value).lower()
                         for c in row)
            if target:
                record["target_rows"].append({"row": index, "cells": {
                    c.coordinate: str(c.value) for c in row if c.value is not None}})
            for cell in row:
                if cell.data_type == "f":
                    record["formula_count"] += 1
                    if len(record["formula_examples"]) < 5 or target:
                        record["formula_examples"].append({"cell": cell.coordinate, "formula": cell.value})
        result["sheets"].append(record)
    workbook.close()
    return result

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source_root", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    sources = list((args.source_root/"研究法/大豆基因/期刊撰寫/table").glob("*FOSCO*.xlsx"))
    sources += list(args.source_root.glob("*.xlsx"))
    sources += list((args.source_root/"研究法/大豆基因/數據/外部資料匯入").glob("*.xlsx"))
    results = []
    for path in sources:
        record = inspect(path)
        results.append(record)
        print(path.name.encode("ascii", errors="backslashreplace").decode(),
              [(s["rows"], s["formula_count"], len(s["target_rows"])) for s in record["sheets"]], flush=True)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(results, ensure_ascii=False, indent=2)+"\n", encoding="utf-8")
