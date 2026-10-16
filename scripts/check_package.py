"""Extract a review ZIP and test its Python and R workflows in a fresh folder.

Prerequisites: Python dependencies installed; Rscript available (or --rscript).
This runs small execution tests and verifies the archived full-run outputs.
It does not repeat the expensive 10,000-draw analyses. Logs are retained.
"""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import time
from zipfile import ZipFile


def main():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", type=Path, default=root.parent/"soybean-cold-tolerance-github.zip")
    parser.add_argument("--rscript", default="Rscript")
    args = parser.parse_args()
    extraction = Path(tempfile.mkdtemp(prefix="dual-language-package-check-", dir=root.parent))
    with ZipFile(args.archive) as archive:
        if archive.testzip() is not None:
            raise ValueError("Archive integrity failed")
        for name in archive.namelist():
            if not (extraction/name).resolve().is_relative_to(extraction.resolve()):
                raise ValueError("Unsafe archive path")
        archive.extractall(extraction)
    workspace = extraction/root.name
    logs = extraction/"logs"
    logs.mkdir()
    py, r = [sys.executable], [args.rscript, "--vanilla"]
    commands = [
        ("r_numerical", r+["tests/test_core.R"]),
        ("r_workflow", r+["tests/test_integration.R"]),
        ("python_tests", py+["-m", "unittest", "discover", "-s", "tests", "-p", "test_*.py"]),
        ("shared_null", py+["tests/check_cross_language.py", "--rscript", args.rscript]),
        ("r_example", r+["examples/walkthrough.R"]),
        ("python_example", py+["examples/walkthrough.py"]),
        ("r_membership", r+["scripts/run_analysis.R"]),
        ("r_plots", r+["scripts/plot_results.R"]),
        ("python_membership", py+["scripts/run_analysis.py"]),
        ("python_plots", py+["scripts/plot_results.py"]),
        ("full_output_verification", py+["scripts/verify_full_runs.py"]),
    ]
    report = {"archive_checkpoint": str(args.archive.resolve()),
              "archive_checkpoint_sha256": hashlib.sha256(args.archive.read_bytes()).hexdigest(),
              "extracted_workspace": str(workspace), "checks": [], "code_and_config_sha256": {}}
    for label, command in commands:
        start = time.perf_counter()
        result = subprocess.run(command, cwd=workspace, capture_output=True,
                                encoding="utf-8", errors="replace")
        (logs/f"{label}.log").write_text(result.stdout+result.stderr, encoding="utf-8")
        report["checks"].append({"name": label, "command": command, "exit_code": result.returncode,
                                 "elapsed_seconds": round(time.perf_counter()-start, 2)})
        print(f"{'PASS' if result.returncode == 0 else 'FAIL'}: {label}", flush=True)
    for pattern in ("R/*.R", "python/*.py", "scripts/*.py", "scripts/*.R", "tests/*.py",
                    "tests/*.R", "examples/*", "config*.R", "config*.json", "requirements.txt"):
        for path in sorted(workspace.glob(pattern)):
            report["code_and_config_sha256"][path.relative_to(workspace).as_posix()] = hashlib.sha256(path.read_bytes()).hexdigest()
    passed = all(item['exit_code'] == 0 for item in report['checks'])
    report['status'] = 'passed' if passed else 'failed'
    output = root/"reports/latest_dual_language_package_check.json"
    output.write_text(json.dumps(report, indent=2)+"\n", encoding="utf-8")
    print(f"Report: {output}\nLogs: {logs}")
    return 0 if passed else 1


if __name__ == "__main__":
    sys.exit(main())
