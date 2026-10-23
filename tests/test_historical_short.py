"""Check recalculation, historical boundaries, and independence from old scores."""
import csv
import hashlib
import math
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"python"))
from ctpathway import read_tsv
from historical_short import (SOURCE_FILES, count_snp_markers, reconstruct_input, resolve_source_files,
    historical_snp_imputation, historical_fosco, historical_score_export, fit_historical_short)


def fixture(folder):
    """Only four upstream sources. No saved inputs or score references exist."""
    paths = {name: folder/relative for name, relative in SOURCE_FILES.items()}
    for path in paths.values():
        path.parent.mkdir(parents=True, exist_ok=True)
    expression = [["index", "index2", "gene", "logFC", "p_short", "FDR", "logFC_middle", "p_middle"]]
    correspondence = [["v1", "v1_1", "v2"]]
    lengths = [["gene", "chr", "start", "end", "length"]]
    blocks = []
    for index in range(20):
        gene = f"Glyma.01g{index+1:06d}"
        length, count = (52247, 0) if index == 19 else (100+index, index % 5+1)
        expression.append([index, index, f"old{index}", .1, .02+index*.035, 1, .2, .8-index*.027])
        correspondence.append([f"old{index}", f"old{index}", gene])
        lengths.append([gene, 1, 1, length+1, length])
        blocks.extend([gene]+[f"AX-{j}" for j in range(count)]+["END", ""])
    for gene, count in [("Glyma.06g180700", 16), ("Glyma.10g123400", 7)]:
        lengths.append([gene, 1, 1, 52211, 52210])
        blocks.extend([gene]+[f"AX-{j}" for j in range(count)]+["END", ""])
    # One unmatched and one ambiguous alias. Matching both aliases in the SAME
    # correspondence row is valid; matching two rows is not.
    expression.extend([[20, 20, "unmapped", .1, .5, 1, .2, .5],
                       [21, 21, "ambiguous", .1, .5, 1, .2, .5]])
    correspondence.extend([["ambiguous", "unused", "Glyma.01g000001"],
                           ["ambiguous", "unused2", "Glyma.01g000002"]])
    for name, rows in [("expression", expression), ("correspondence", correspondence), ("gene_length", lengths)]:
        with paths[name].open("w", encoding="utf-8", newline="") as handle:
            csv.writer(handle, delimiter="\t").writerows(rows)
    paths["snp_markers"].write_text("\n".join(blocks), encoding="utf-8")
    return paths


class HistoricalShortTests(unittest.TestCase):
    def test_portable_manifest_and_corruption_do_not_fall_back(self):
        """Moving the input files must work; corruption must fail in both languages."""
        with tempfile.TemporaryDirectory() as temporary:
            folder = Path(temporary)
            paths = fixture(folder)
            expected, _, _ = reconstruct_input(folder)
            portable = folder / "portable inputs with spaces"
            portable.mkdir()
            entries = []
            for key, source in paths.items():
                target = portable / (key + ".txt")
                shutil.copyfile(source, target)
                contents = target.read_bytes()
                entries.append(dict(source=key, file=target.name,
                                    sha256=hashlib.sha256(contents).hexdigest(),
                                    md5=hashlib.md5(contents).hexdigest()))
            with (portable / "inputs.tsv").open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, delimiter="\t", fieldnames=list(entries[0]))
                writer.writeheader()
                writer.writerows(entries)
            self.assertEqual(reconstruct_input(portable)[0], expected)
            rscript = shutil.which("Rscript") or "C:/Program Files/R/R-4.5.1/bin/Rscript.exe"
            rcheck = [rscript, "--vanilla", "-e",
                      'source("R/historical_short.R"); x <- reconstruct_historical_input(commandArgs(TRUE)[1]); stopifnot(nrow(x$rows) == 20L)',
                      str(portable)]
            if Path(rscript).is_file():
                result = subprocess.run(rcheck, cwd=ROOT, capture_output=True, text=True)
                self.assertEqual(result.returncode, 0, result.stderr)
            with (portable / "expression.txt").open("a", encoding="utf-8") as handle:
                handle.write("corruption\n")
            with self.assertRaisesRegex(ValueError, "checksum mismatch"):
                resolve_source_files(portable)
            if Path(rscript).is_file():
                result = subprocess.run(rcheck, cwd=ROOT, capture_output=True, text=True)
                self.assertNotEqual(result.returncode, 0)
                self.assertIn("checksum mismatch", result.stderr)

    def test_original_bins_include_zero_and_leave_no_donor_bin_unresolved(self):
        rows = [{"source_row": i+2, "gene_id": f"g{i}", "gene_length": length, "snp_markers": snp}
                for i, (length, snp) in enumerate([(1, 4), (2000, 0), (2001, 0), (98001, 0)])]
        values, audit = historical_snp_imputation(rows)
        np.testing.assert_array_equal(values, [4, 2, 0, 0])
        self.assertEqual(audit[0]["mean_including_zeros"], 2)

    def test_export_precision_zero_floor_and_rounding_boundary(self):
        self.assertEqual(historical_score_export(1)[2], "1.00E-16")
        self.assertEqual(historical_score_export(0)[2], "Inf")
        self.assertEqual(historical_score_export(.9947477)[2], "0.002287056")
        self.assertEqual(historical_score_export(1.176836406102666e-14)[2], "13.9292839")

    def test_formula_and_input_validation(self):
        self.assertAlmostEqual(historical_fosco([.2], [4], .5)[0], .36)
        self.assertEqual(historical_fosco([.2], [0], .5)[0], 0)
        self.assertNotEqual(historical_fosco([1.14e-14], [10], .01)[0],
                            -np.expm1(10**.01*np.log1p(-1.14e-14)))
        for p, n in [([0, .3], [1, 2]), ([.2, .3], [1, -1]), ([.2, .3], [1, 1.5])]:
            with self.assertRaises(ValueError):
                fit_historical_short(p, n)

    def check_cli(self, language, executable):
        scratch = ROOT/"tests/tmp"
        scratch.mkdir(exist_ok=True)
        with tempfile.TemporaryDirectory(dir=scratch) as temp:
            folder = Path(temp)
            paths = fixture(folder)
            rows, exclusions, _ = reconstruct_input(folder)
            self.assertEqual(len(rows), 20)
            self.assertEqual({r["reason"] for r in exclusions}, {"no_correspondence", "ambiguous_correspondence_rows"})
            script = ROOT/f"scripts/reproduce_historical_short.{language}"
            command = [executable, str(script), str(folder)]
            strict = subprocess.run(command+[str(folder/"strict")], cwd=ROOT, capture_output=True, text=True)
            self.assertEqual(strict.returncode, 2, strict.stdout+strict.stderr)
            self.assertEqual(len(read_tsv(folder/"strict/unresolved_rows.tsv")), 1)

            def run(label):
                result = subprocess.run(command+[str(folder/label), "--zero-repair-hypothesis"],
                                        cwd=ROOT, capture_output=True, text=True)
                self.assertEqual(result.returncode, 0, result.stdout+result.stderr)
                return read_tsv(folder/label/"recalculated_short.tsv")
            baseline = run("hypothesis")
            repaired = baseline[-1]
            self.assertEqual(float(repaired["snp"]), 12)
            audit = read_tsv(folder/"hypothesis/zero_repair_hypothesis.tsv")[0]
            self.assertEqual(float(audit["donor_mean"]), 11.5)
            self.assertEqual(audit["evidence_status"], "unconfirmed_reconstruction_hypothesis")
            # Any accidental attempt to read a saved probability/score fails.
            for name in ["gene_p_len.txt", "Yamasaki_Go校正後1213.txt", "Yamasaki_Go校正後1213_G_log10.txt"]:
                (paths["snp_markers"].parent/name).write_text("THIS IS NOT A SCORE TABLE", encoding="utf-8")
            self.assertEqual(baseline, run("changed_references"))
            with paths["expression"].open(encoding="utf-8", newline="") as handle:
                source = list(csv.reader(handle, delimiter="\t"))
            source[20][4] = ".1"  # Change only the exceptional row's raw short input.
            with paths["expression"].open("w", encoding="utf-8", newline="") as handle:
                csv.writer(handle, delimiter="\t").writerows(source)
            changed = run("changed_input")
            self.assertNotEqual(repaired["score_short"], changed[-1]["score_short"])
            self.assertEqual(float(repaired["snp"]), float(changed[-1]["snp"]))

    def test_python_cli_recalculates_without_reference_dependency(self):
        self.check_cli("py", sys.executable)

    def test_r_cli_recalculates_without_reference_dependency(self):
        rscript = shutil.which("Rscript") or "C:/Program Files/R/R-4.5.1/bin/Rscript.exe"
        if not Path(rscript).is_file():
            self.skipTest("Rscript unavailable; standalone R unit tests run in the R CI job")
        self.check_cli("R", rscript)


if __name__ == "__main__":
    unittest.main()
