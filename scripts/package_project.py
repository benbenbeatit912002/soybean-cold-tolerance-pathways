"""Build a review ZIP including local results, excluding temporary test files."""
from pathlib import Path
from zipfile import ZipFile, ZIP_DEFLATED

root = Path(__file__).resolve().parents[1]
archive = root.parent / "soybean-cold-tolerance-github.zip"
included = []
with ZipFile(archive, "w", compression=ZIP_DEFLATED, compresslevel=6) as z:
    for path in sorted(root.rglob("*")):
        relative = path.relative_to(root)
        if not path.is_file() or any(p in (".git", "__pycache__") or p.startswith(".venv") for p in relative.parts):
            continue
        if relative.parts[:2] == ("tests", "tmp") or path.suffix == ".pyc":
            continue
        z.write(path, (Path(root.name) / relative).as_posix())
        included.append(relative.as_posix())
with ZipFile(archive) as z:
    assert z.testzip() is None
    assert len(z.namelist()) == len(included)
print(f"Created {archive.name}: {len(included)} files, {archive.stat().st_size:,} bytes; archive integrity passed.")
