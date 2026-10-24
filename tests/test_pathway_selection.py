"""Selection replay is checked against references only after calculation."""
import copy
import csv
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT / "python"))
sys.path.insert(0,str(ROOT / "scripts"))
from pathway_selection import select_candidates
from verify_presoynet import compare


class CandidateSelectionTests(unittest.TestCase):
    def setUp(self):
        self.config = json.loads((ROOT / "config.selection-candidate.json").read_text())
        with (ROOT / "data/reference/crosstalk/selection_candidates.tsv").open(encoding="utf-8",newline="") as f:
            self.rows = list(csv.DictReader(f,delimiter="\t"))

    def test_reconstruction_without_target_ids(self):
        audit=select_candidates(self.rows,self.config)
        with (ROOT / "data/reference/crosstalk/selected_pathways.tsv").open(encoding="utf-8",newline="") as f:
            expected={(r["phase"],r["go_id"]) for r in csv.DictReader(f,delimiter="\t")}
        self.assertEqual({(r["phase"],r["go_id"]) for r in audit if r["selected"]},expected)
        # Changing a scientifically meaningful boundary must change the result.
        inclusive=copy.deepcopy(self.config); inclusive["phases"]["short"]["min_ct"]=5
        changed=select_candidates(self.rows,inclusive)
        self.assertEqual(sum(r["selected"] for r in changed if r["phase"]=="short"),16)
        bh=copy.deepcopy(self.config); bh["probability_field"]="p_bh"
        result=select_candidates(self.rows,bh)
        self.assertFalse(next(r["selected"] for r in result if r["phase"]=="short" and r["go_id"]=="GO:0009737"))
        with self.assertRaises(ValueError): select_candidates(self.rows+[self.rows[0]],self.config)

    def test_python_r_relocated_cli(self):
        rscript=shutil.which("Rscript") or "C:/Program Files/R/R-4.5.1/bin/Rscript.exe"
        if not Path(rscript).is_file(): self.skipTest("Rscript is not installed")
        with tempfile.TemporaryDirectory(prefix="selection with spaces ") as folder:
            destination=Path(folder)
            for language,exe,extension in (("python",[sys.executable],"py"),("r",[rscript,"--vanilla"],"R")):
                config=ROOT / ("config.selection-candidate.json" if language=="python" else "config.selection-candidate.R")
                command=[*exe,str(ROOT/f"scripts/select_pathways.{extension}"),str(ROOT/"data/reference/crosstalk/selection_candidates.tsv"),str(config),str(destination/language)]
                run=subprocess.run(command,cwd=destination,capture_output=True,text=True)
                self.assertEqual(run.returncode,0,run.stdout+run.stderr)
            compare(destination/"python/selection_audit.tsv",destination/"r/selection_audit.tsv")
            compare(destination/"python/selected_pathways.tsv",destination/"r/selected_pathways.tsv")


if __name__=="__main__": unittest.main()
