"""Build a portable source ZIP; generated results are optional release assets.

Default: python scripts/package_project.py
Optional result bundle: python scripts/package_project.py --include-results
No Git checkout, environment, cache, or test scratch directory is packaged.
"""
import argparse
import hashlib
from pathlib import Path
from zipfile import ZipFile, ZIP_DEFLATED

ROOT = Path(__file__).resolve().parents[1]
ARCHIVE_ROOT = "soybean-cold-tolerance-pathways"
FOLDERS = (".github", "R", "python", "scripts", "tests", "examples", "data", "docs", "reports", "legacy")


def package_files(include_results=False):
    for path in sorted(ROOT.iterdir()):
        if path.is_file() and (path.suffix in (".py", ".R", ".json", ".md", ".cff") or
                               path.name in ("requirements.txt", ".gitignore", ".gitattributes", "LICENSE")):
            yield path
    for folder in (*FOLDERS, *(("results",) if include_results else ())):
        for path in sorted((ROOT / folder).rglob("*")):
            relative = path.relative_to(ROOT)
            if (not path.is_file() or path.suffix == ".pyc" or "__pycache__" in relative.parts or
                    relative.parts[:2] == ("tests", "tmp") or relative.parts[:2] == ("results", "package_checks")):
                continue
            yield path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--include-results", action="store_true")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    kind = "with-results" if args.include_results else "source"
    output = args.output or ROOT / "dist" / f"{ARCHIVE_ROOT}-{kind}.zip"
    output.parent.mkdir(parents=True, exist_ok=True)
    paths = list(package_files(args.include_results))
    with ZipFile(output, "w", compression=ZIP_DEFLATED, compresslevel=6) as archive:
        for path in paths:
            archive.write(path, (Path(ARCHIVE_ROOT) / path.relative_to(ROOT)).as_posix())
    with ZipFile(output) as archive:
        if archive.testzip() is not None or len(archive.namelist()) != len(paths):
            raise ValueError("Archive integrity failed")
    digest = hashlib.sha256(output.read_bytes()).hexdigest()
    output.with_suffix(output.suffix + ".sha256").write_text(f"{digest}  {output.name}\n", encoding="utf-8")
    print(f"Created {output}: {len(paths)} files, {output.stat().st_size:,} bytes; integrity passed.")


if __name__ == "__main__":
    main()
