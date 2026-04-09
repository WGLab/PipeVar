#!/usr/bin/env python3
import argparse
import csv
from pathlib import Path


def count_rows(path: Path) -> int:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return max(sum(1 for _ in handle) - 1, 0)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--raw-dir", required=True)
    parser.add_argument("--prepared-dir", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    raw_dir = Path(args.raw_dir)
    prepared_dir = Path(args.prepared_dir)
    out_path = Path(args.out)

    rows = []
    for root, label in ((raw_dir, "raw"), (prepared_dir, "prepared")):
        for path in sorted(root.glob("*")):
            if path.is_file():
                rows.append(
                    {
                        "layer": label,
                        "filename": path.name,
                        "bytes": str(path.stat().st_size),
                        "rows": str(count_rows(path)) if path.suffix in {".tsv", ".txt"} else ".",
                    }
                )

    with out_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["layer", "filename", "bytes", "rows"], delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
