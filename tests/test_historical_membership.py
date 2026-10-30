"""Validate recalculated membership independently of saved pathway statistics."""
import csv
from decimal import Decimal
import math
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
REFERENCE = ROOT / "data/reference/historical_membership"


def read(path):
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def printed_tolerance(value):
    if float(value) in (0, 1):
        return 0.0 if float(value) == 0 else 1e-14
    return float(Decimal(10) ** Decimal(value).as_tuple().exponent) / 2 + 1e-18


class HistoricalMembershipTests(unittest.TestCase):
    def test_independent_recount_and_saved_precision(self):
        rscript = shutil.which("Rscript") or "C:/Program Files/R/R-4.5.1/bin/Rscript.exe"
        languages = [("python", [sys.executable], "py")]
        if Path(rscript).is_file():
            languages.append(("r", [rscript, "--vanilla"], "R"))
        with tempfile.TemporaryDirectory(prefix="membership audit with spaces ") as folder:
            output = Path(folder)
            results = {}
            for language, executable, extension in languages:
                command = executable + [str(ROOT / f"scripts/audit_historical_membership.{extension}"),
                    str(ROOT / "data/reference/historical_input.tsv"), str(output / language)]
                run = subprocess.run(command, cwd=output, capture_output=True, text=True)
                self.assertEqual(run.returncode, 0, run.stdout + run.stderr)
                for phase in ("short", "middle"):
                    membership = {r["go_id"]: r for r in read(output / language / f"{phase}_membership.tsv")}
                    for row in membership.values():
                        excluded = row["excluded_ct_ids"].split(";") if row["excluded_ct_ids"] else []
                        self.assertEqual(int(row["primary_ct_count"])-int(row["measured_ct_count"]), len(excluded))
                    actual = {r["go_id"]: r for r in read(output / language / f"{phase}_historical_hypergeometric.tsv")}
                    saved = {r["go_id"]: r for r in read(REFERENCE / f"{phase}_0109.tsv")}
                    self.assertEqual(set(actual), set(saved))
                    for go, row in saved.items():
                        self.assertEqual(int(actual[go]["measured_ct_count"]), int(row["ct_count"]))
                        for old_field, new_field in (("p_value", "one_minus_cdf"), ("p_bh", "one_minus_cdf_bh")):
                            # R reproduces the archived subtraction arithmetic.
                            # Python can differ by one CDF step near machine epsilon.
                            tolerance = printed_tolerance(row[old_field])
                            if language == "python" and phase == "middle":
                                tolerance += 2e-16 if old_field == "p_value" else 1e-15
                            self.assertLessEqual(abs(float(actual[go][new_field])-float(row[old_field])), tolerance)
                    results[language, phase] = actual
                missing = {r["gene_id"] for r in read(output / language / "missing_ctgenes.tsv") if r["phase"] == "short"}
                primary = {r["gene_id"] for r in read(ROOT / "data/input/ct_short.tsv")}
                historical_scored = {r["gene_id"] for r in read(REFERENCE / "short_scored_genes.tsv")}
                self.assertEqual(primary-missing, historical_scored)
            if len(languages) == 2:
                for phase in ("short", "middle"):
                    for go, row in results["python", phase].items():
                        for field in ("p_value", "p_bh"):
                            self.assertTrue(math.isclose(float(row[field]), float(results["r", phase][go][field]),
                                                         rel_tol=2e-11, abs_tol=0))


if __name__ == "__main__":
    unittest.main()
