"""Validate reference matching for independent audits; reject obsolete score copying."""
from copy import deepcopy
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python"))
from ctpathway import (fit_fosco, impute_snp_counts, match_short_reference,
                       read_tsv, write_tsv)


def fixture():
    """Keep duplicate gene IDs distinct through their original source-row keys."""
    rows = [{"source_row": str(i + 2),
             "gene_id": f"Glyma.01G{max(i, 1):06d}",
             "raw_short": str(.02 + i * .035),
             "raw_middle": str(.8 - i * .027),
             "gene_length": str(52247 if i == 19 else 100 + i),
             "snp": str(0 if i == 19 else i % 5 + 1)}
            for i in range(20)]
    references = [{**row, "gene_id": row["gene_id"].replace("Glyma.01G", "Glyma.01g"),
                   "reference_adjusted_short": "0.9000000000000000",
                   # Extra stored digits must not be lost to a float round-trip.
                   "reference_score_short": "0.045757490560675125409944193489",
                   "reference_snp_short": "12",
                   "reference_log_row": str(100 + i)}
                  for i, row in enumerate(rows)]
    return rows, references


class ShortReferenceTest(unittest.TestCase):
    def test_matches_shuffled_rows_and_preserves_stored_text(self):
        rows, reference = fixture()
        original = deepcopy(rows)
        matched = match_short_reference(rows, list(reversed(reference)))
        self.assertEqual([r["source_row"] for r in matched],
                         [r["source_row"] for r in rows])
        self.assertEqual([r["reference_log_row"] for r in matched],
                         [str(100 + i) for i in range(20)])
        self.assertEqual(len(matched), 20)
        self.assertEqual(rows[0]["gene_id"], rows[1]["gene_id"])
        self.assertNotEqual(matched[0]["reference_log_row"],
                            matched[1]["reference_log_row"])
        self.assertTrue(all(r["reference_score_short"] == reference[0]["reference_score_short"]
                            for r in matched))
        self.assertEqual(rows, original, "Matching must not modify the correction input")

    def test_rejects_mismatched_membership_and_identity(self):
        rows, reference = fixture()
        cases = []
        cases.append((rows, reference[:-1]))
        cases.append((rows, reference + [dict(reference[-1], source_row="999")]))
        for field, bad in [("source_row", "999"), ("source_row", "2"),
                           ("source_row", "21.5"), ("source_row", "0"),
                           ("source_row", "nan"), ("gene_id", "Glyma.99G999999"),
                           ("raw_short", "0.999"), ("gene_length", "52248"),
                           ("snp", "1")]:
            altered = deepcopy(reference)
            altered[-1][field] = bad
            cases.append((rows, altered))
        for field, bad in [("source_row", "2"), ("source_row", "21.5"),
                           ("source_row", "0"), ("source_row", "inf")]:
            altered = deepcopy(rows)
            altered[-1][field] = bad
            cases.append((altered, reference))
        missing = deepcopy(reference)
        del missing[0]["reference_score_short"]
        cases.append((rows, missing))
        missing = deepcopy(rows)
        del missing[0]["source_row"]
        cases.append((missing, reference))
        for index, (left, right) in enumerate(cases):
            with self.subTest(case=index):
                with self.assertRaises(ValueError):
                    match_short_reference(left, right)

    def test_rejects_invalid_reference_values(self):
        rows, reference = fixture()
        invalid = {
            "reference_adjusted_short": ["0", "-0.1", "1.1", "nan", "inf", ""],
            "reference_score_short": ["-0.1", "nan", "inf", ""],
            "reference_snp_short": ["0", "-1", "2.5", "nan", "inf", ""],
            "reference_log_row": ["0", "-1", "2.5", "nan", "inf", ""],
        }
        for field, values in invalid.items():
            for value in values:
                with self.subTest(field=field, value=value):
                    altered = deepcopy(reference)
                    altered[0][field] = value
                    with self.assertRaises(ValueError):
                        match_short_reference(rows, altered)

    def test_cli_rejects_obsolete_score_copy_option(self):
        result = subprocess.run([sys.executable, str(ROOT / "scripts/correct_scores.py"),
                                 "unused.tsv", "unused-output", "--short-reference", "unused-reference.tsv"],
                                cwd=ROOT, capture_output=True, text=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("unrecognized arguments", result.stderr)


if __name__ == "__main__":
    unittest.main()
