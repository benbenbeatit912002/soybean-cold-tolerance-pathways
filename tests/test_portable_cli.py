"""Check the documented demo from another directory and protect existing results."""
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]


class PortableEntryPoints(unittest.TestCase):
    def check_entry(self, prefix, suffix):
        with tempfile.TemporaryDirectory(prefix="caller with spaces ") as temporary:
            caller = Path(temporary)
            command = [*prefix, str(ROOT / f"run_pipeline.{suffix}"), "--mode", "demo",
                       "--output", "relative results"]
            completed = subprocess.run(command, cwd=caller, capture_output=True, text=True)
            self.assertEqual(completed.returncode, 0, completed.stdout + completed.stderr)
            output = caller / "relative results/analysis"
            self.assertEqual(len(list(output.glob("*.png"))), 6)
            self.assertFalse((output / "short_maxmean.tsv").exists())
            table = output / "short_sumstat.tsv"
            original = table.read_bytes()
            again = subprocess.run(command, cwd=caller, capture_output=True, text=True)
            self.assertNotEqual(again.returncode, 0)
            self.assertIn("new or empty", again.stdout + again.stderr)
            self.assertEqual(table.read_bytes(), original, "A repeat invocation must not alter existing results")

    def test_python_from_another_directory(self):
        self.check_entry([sys.executable], "py")

    def test_r_from_another_directory(self):
        rscript = shutil.which("Rscript") or "C:/Program Files/R/R-4.5.1/bin/Rscript.exe"
        if not Path(rscript).is_file():
            self.skipTest("Rscript is not installed")
        self.check_entry([rscript, "--vanilla"], "R")


if __name__ == "__main__":
    unittest.main()
