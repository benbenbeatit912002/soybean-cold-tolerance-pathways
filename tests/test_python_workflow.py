"""Run the public Python analysis, correction, and plotting commands on toy data."""
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"python"))
from ctpathway import read_tsv, write_tsv


class WorkflowTest(unittest.TestCase):
    def command(self, script, *arguments):
        result = subprocess.run([sys.executable, str(ROOT/script), *map(str, arguments)],
                                cwd=ROOT, capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stdout+result.stderr)

    def test_both_phases_all_methods_correction_and_plots(self):
        scratch = ROOT/"tests/tmp"
        scratch.mkdir(exist_ok=True)
        with tempfile.TemporaryDirectory(prefix="python-workflow-", dir=scratch) as temporary:
            folder = Path(temporary)
            inputs, output = folder/"input", folder/"results"
            inputs.mkdir()
            genes = [f"Glyma.01G{i:06d}" for i in range(1, 21)]
            for phase in ("short", "middle"):
                write_tsv([{"gene_id": gene, "prior_score": 8} for gene in genes[:4]], inputs/f"ct_{phase}.tsv")
                write_tsv([dict(gene_id=gene, **{f"source_col_{c}": i for c in range(2, 7)})
                           for i, gene in enumerate(genes, start=1)], inputs/f"external_{phase}.tsv")
            (inputs/"pathways.gmt").write_text("\n".join("\t".join([go, "toy pathway", *members])
                for go, members in [("GO:1", genes[:6]), ("GO:2", genes[-6:])])+"\n")
            config = json.loads((ROOT/"config.python.json").read_text())
            config.update(input_dir=str(inputs), output_dir=str(output), score_column=6,
                          score_kind="score", permutations=39, methods=["sumstat", "gsea", "maxmean"])
            config_path = folder/"config.json"
            config_path.write_text(json.dumps(config))
            self.command("scripts/run_analysis.py", config_path)
            for phase in config['phases']:
                for method in config['methods']:
                    rows = read_tsv(output/f"{phase}_{method}.tsv")
                    self.assertEqual(len(rows), 2)
                    self.assertTrue(all(float(r['p_value']) >= 1/40 for r in rows))
            self.command("scripts/plot_results.py", output)
            self.assertEqual(len(list(output.glob("*.png"))), 8)
            raw = folder/"raw.tsv"
            write_tsv([{"gene_id": gene, "raw_short": .01+i*.04,
                        "raw_middle": .8-i*.03, "snp": i % 5+1}
                       for i, gene in enumerate(genes)], raw)
            correction = folder/"correction"
            self.command("scripts/correct_scores.py", raw, correction)
            corrected = read_tsv(correction/"corrected_scores.tsv")
            self.assertEqual(len(corrected), 20)
            self.assertTrue(all(0 < float(row[f'adjusted_{phase}']) <= 1
                                for row in corrected for phase in ('short', 'middle')))
            self.assertEqual(len(list(correction.glob("*.png"))), 2)


if __name__ == "__main__":
    unittest.main()
