"""Extract a source ZIP and exercise both languages without the original folder.

Default: unit/integration tests, both demos and both real-data quick runs.
--full: also rerun 10,000-draw analyses, including experimental MaxMean.
The clean extraction, command logs and exit codes are retained under results/.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import time
from zipfile import ZipFile

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", type=Path, default=ROOT / "dist/soybean-cold-tolerance-pathways-source.zip")
    parser.add_argument("--rscript", default="Rscript")
    parser.add_argument("--full", action="store_true")
    args = parser.parse_args()
    scratch = ROOT / "results/package_checks"
    scratch.mkdir(parents=True, exist_ok=True)
    extraction = Path(tempfile.mkdtemp(prefix="download test ", dir=scratch))
    with ZipFile(args.archive) as archive:
        if archive.testzip() is not None:
            raise ValueError("Archive integrity failed")
        for name in archive.namelist():
            if not (extraction / name).resolve().is_relative_to(extraction.resolve()):
                raise ValueError("Unsafe archive path")
        archive.extractall(extraction)
    workspace = extraction / "soybean-cold-tolerance-pathways"
    if (workspace / "results").exists():
        raise ValueError("Use the source-only ZIP: preexisting outputs could hide missing calculation steps")
    logs = extraction / "logs"
    logs.mkdir()
    py, r = [sys.executable], [args.rscript, "--vanilla"]
    report = dict(archive=str(args.archive.resolve()), archive_sha256=hashlib.sha256(args.archive.read_bytes()).hexdigest(),
                  extracted_workspace=str(workspace), checks=[], full_simulations_requested=args.full,
                  original_project_used=False, status="running", tested_file_sha256={})
    for path in workspace.rglob("*"):
        if path.is_file():
            report["tested_file_sha256"][path.relative_to(workspace).as_posix()] = hashlib.sha256(path.read_bytes()).hexdigest()
    report_path = ROOT / "reports/latest_portable_package_check.json"

    def save():
        report_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")

    def execute(label, command, expected=0, outside=False):
        start = time.perf_counter()
        with (logs / f"{label}.log").open("w", encoding="utf-8") as handle:
            result = subprocess.run(command, cwd=extraction if outside else workspace, stdout=handle, stderr=subprocess.STDOUT)
        item = dict(name=label, command=command, exit_code=result.returncode, expected_exit_code=expected,
                    elapsed_seconds=round(time.perf_counter() - start, 2), log=str(logs / f"{label}.log"))
        print(f"{'PASS' if result.returncode == expected else 'FAIL'}: {label} ({item['elapsed_seconds']} s)", flush=True)
        return item

    commands = [("python_tests", py + ["-m", "unittest", "discover", "-s", "tests", "-p", "test_*.py"]),
                *[(name, r + [f"tests/{name}.R"]) for name in ("test_core", "test_integration", "test_short_reference", "test_historical_short")],
                ("shared_null", py + ["tests/check_cross_language.py", "--rscript", args.rscript])]
    save()
    for label, command in commands:
        report["checks"].append(execute(label, command))
        save()
    # R and Python have independent output folders. Every entry-point test starts
    # outside the repository, using a relocated path that contains spaces.
    with ThreadPoolExecutor(max_workers=2) as pool:
        for mode in ("demo", "quick", *(("full",) if args.full else ())):
            jobs = []
            for language, executable, extension in (("python", py, "py"), ("r", r, "R")):
                command = executable + [str(workspace / f"run_pipeline.{extension}"), "--mode", mode, "--include-maxmean"]
                if mode != "demo":
                    command += ["--use-reconstruction-hypotheses"]
                jobs.append(pool.submit(execute, f"{language}_{mode}", command, outside=True))
            report["checks"].extend(job.result() for job in jobs)
            save()
    # Deliberately omit the hypotheses: unresolved scores must stop before enrichment.
    for language, executable, extension in (("python", py, "py"), ("r", r, "R")):
        command = executable + [str(workspace / f"run_pipeline.{extension}"), "--mode", "quick", "--output", str(workspace / f"results/{language}_strict")]
        report["checks"].append(execute(f"{language}_strict", command, expected=2, outside=True))
    mode = "full" if args.full else "quick"
    report["checks"].append(execute("portable_results", py + ["scripts/verify_portable_runs.py",
        f"results/python_{mode}", f"results/r_{mode}", "--permutations", "10000" if args.full else "199"]))
    report["status"] = "passed" if all(c["exit_code"] == c["expected_exit_code"] for c in report["checks"]) else "failed"
    save()
    print(f"Report: {report_path}\nLogs: {logs}", flush=True)
    return 0 if report["status"] == "passed" else 1


if __name__ == "__main__":
    sys.exit(main())
