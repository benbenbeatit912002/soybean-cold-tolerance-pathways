"""Maintainer utility: copy four original inputs into the portable data layout.

End users already have these files and do not need the original project folder.
Copies retain every byte. The existing correspondence copy is reused only when
its checksum agrees. Original files are never modified.
"""
import argparse
import csv
import hashlib
from pathlib import Path
import shutil
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python"))
from historical_short import SOURCE_FILES


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("original_project", type=Path)
    args = parser.parse_args()
    destination = ROOT / "data/upstream"
    destination.mkdir(parents=True, exist_ok=True)
    names = {"expression": "expression.tsv", "correspondence": "../input/gene_correspondence.tsv",
             "gene_length": "gene_lengths.tsv", "snp_markers": "snp_markers.txt"}
    entries = []
    for name, relative in SOURCE_FILES.items():
        source = args.original_project / relative
        contents = source.read_bytes()
        target = destination / names[name]
        if target.exists() and target.read_bytes() != contents:
            raise ValueError(f"Refusing to replace a different staged input: {target}")
        if not target.exists():
            shutil.copyfile(source, target)
        entries.append(dict(source=name, file=names[name], bytes=len(contents),
                            sha256=hashlib.sha256(contents).hexdigest(), md5=hashlib.md5(contents).hexdigest(),
                            original_relative_path=relative, transformation="Byte-identical copy"))
    with (destination / "inputs.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(entries[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(entries)
    print(f"Staged {len(entries)} upstream sources; hashes saved to data/upstream/inputs.tsv")


if __name__ == "__main__":
    main()
