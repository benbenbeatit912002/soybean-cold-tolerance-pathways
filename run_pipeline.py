"""Run a complete Python workflow from any working directory.

Examples (from the repository folder):
    python run_pipeline.py --mode demo
    python run_pipeline.py --mode full --use-reconstruction-hypotheses

Demo uses invented data. Quick/full rebuild real scores from included upstream
inputs. Reconstruction hypotheses are explicit and remain historically unproven.
"""
import argparse
import json
from pathlib import Path
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parent


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", choices=("demo", "quick", "full"), default="demo")
    parser.add_argument("--output", type=Path, help="New output directory; relative to your current directory")
    parser.add_argument("--use-reconstruction-hypotheses", action="store_true",
                        help="Enable the two documented, unconfirmed historical score repairs")
    parser.add_argument("--include-maxmean", action="store_true", help="Also run experimental MaxMean")
    args = parser.parse_args()
    output = (args.output or ROOT / f"results/python_{args.mode}").resolve()
    if output.exists() and (not output.is_dir() or any(output.iterdir())):
        parser.error("Output must be new or empty. Use --output with a different directory.")
    if args.mode == "demo" and args.use_reconstruction_hypotheses:
        parser.error("Demo uses invented scores and does not accept reconstruction hypotheses.")
    output.mkdir(parents=True, exist_ok=True)
    metadata = dict(language="Python", mode=args.mode, status="running", steps=[],
                    hypotheses_enabled=args.use_reconstruction_hypotheses,
                    historical_repairs_confirmed=False, paper_reproduction_verified=False)

    def save():
        (output / "run_manifest.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")

    def run(script, *arguments):
        command = [sys.executable, str(ROOT / script), *map(str, arguments)]
        print(f"\nRunning {script}", flush=True)
        started = time.perf_counter()
        result = subprocess.run(command, cwd=ROOT)
        metadata["steps"].append(dict(script=script, command=command, exit_code=result.returncode,
                                      elapsed_seconds=round(time.perf_counter() - started, 3)))
        save()
        if result.returncode:
            metadata["status"] = "unresolved_scores" if result.returncode == 2 and "recalculate" in script else "failed"
            save()
            print(f"Stopped. Inspect {output}; no downstream success is being reported.", file=sys.stderr)
            raise SystemExit(result.returncode)

    # 1. Rebuild scores using only packaged raw inputs; demo scores are synthetic.
    save()
    if args.mode == "demo":
        inputs = ROOT / "examples/data"
        scores = inputs / "scores.tsv"
    else:
        inputs = ROOT / "data/input"
        flags = ["--use-reconstruction-hypotheses"] if args.use_reconstruction_hypotheses else []
        run("scripts/recalculate_phase_scores.py", ROOT / "data/upstream", output / "scores", *flags)
        scores = output / "scores/recalculated_scores.tsv"

    # 2. Use named columns: no inference about the older WITHOUTINF score columns.
    config = dict(input_dir=str(inputs), score_input=str(scores), score_field="score_{phase}",
                  score_column=None, score_kind="score", phases=["short", "middle"],
                  min_size=2 if args.mode == "demo" else 6,
                  max_size=18 if args.mode == "demo" else 1500,
                  permutations=10000 if args.mode == "full" else 199, seed=20261021,
                  duplicate_policy="mean", methods=["sumstat", "gsea"],
                  output_dir=str(output / "analysis"),
                  score_provenance="Invented teaching data" if args.mode == "demo" else
                                   "Recalculated historical p-value branch with explicit repair hypotheses",
                  paper_reproduction_verified=False)
    if args.include_maxmean:
        config["methods"].append("maxmean")
    config_path = output / "analysis_config.json"
    config_path.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")

    # 3. Run overlap and score tests, then plot the resulting tables.
    run("scripts/run_analysis.py", config_path)
    run("scripts/plot_results.py", output / "analysis")
    metadata.update(status="completed", permutations=config["permutations"], methods=["hypergeometric", *config["methods"]])
    save()
    print(f"\nCompleted {args.mode}. Tables and plots: {output / 'analysis'}", flush=True)
    if args.mode != "full":
        print("This is an execution/teaching run; use full mode for the 10,000-draw analysis.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
