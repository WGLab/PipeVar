#!/usr/bin/env python3
"""Add sex from a separate metadata CSV to a PipeVar manifest."""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path
from typing import Iterable, Sequence


class EnrichmentError(ValueError):
    """Raised when sex metadata cannot be added unambiguously."""


def positive_column(value: str) -> int:
    try:
        column = int(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("column must be a positive integer") from exc
    if column < 1:
        raise argparse.ArgumentTypeError("column must be a positive integer")
    return column


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Match a PipeVar manifest to a separate metadata CSV and add a "
            "normalized sex column. Column numbers are one-based."
        )
    )
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--metadata", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--manifest-sample-column", type=positive_column, default=1)
    parser.add_argument("--metadata-sample-column", type=positive_column, default=2)
    parser.add_argument("--metadata-sex-column", type=positive_column, default=10)
    parser.add_argument("--strip-suffix", default=".norm.alt")
    return parser.parse_args(argv)


def read_csv(path: Path, label: str) -> tuple[list[str], list[tuple[int, list[str]]]]:
    if not path.is_file():
        raise EnrichmentError(f"{label} CSV not found: {path}")

    try:
        with path.open("r", encoding="utf-8-sig", newline="") as handle:
            reader = csv.reader(handle, strict=True)
            records = [
                (line_number, row)
                for line_number, row in enumerate(reader, start=1)
                if row
            ]
    except (OSError, UnicodeError, csv.Error) as exc:
        raise EnrichmentError(f"could not parse {label} CSV '{path}': {exc}") from exc

    if not records:
        raise EnrichmentError(f"{label} CSV is empty: {path}")

    header = records[0][1]
    expected_columns = len(header)
    rows: list[tuple[int, list[str]]] = []
    for line_number, row in records[1:]:
        if len(row) != expected_columns:
            raise EnrichmentError(
                f"{label} CSV row {line_number} has {len(row)} columns; "
                f"expected {expected_columns}"
            )
        rows.append((line_number, row))
    return header, rows


def column_index(header: Sequence[str], one_based_column: int, label: str) -> int:
    index = one_based_column - 1
    if index >= len(header):
        raise EnrichmentError(
            f"{label} column {one_based_column} is outside the "
            f"{len(header)}-column CSV"
        )
    return index


def normalize_manifest_sample(raw_sample: str, suffix: str) -> str:
    sample = raw_sample.strip()
    if suffix and sample.endswith(suffix):
        sample = sample[: -len(suffix)]
    return sample.strip()


def normalize_sex(raw_sex: str, source: str) -> str:
    sex = raw_sex.strip().lower()
    if sex not in {"male", "female"}:
        raise EnrichmentError(
            f"invalid sex '{raw_sex}' for {source}; expected MALE or FEMALE"
        )
    return sex


def build_sex_lookup(
    header: Sequence[str],
    rows: Iterable[tuple[int, list[str]]],
    sample_column: int,
    sex_column: int,
) -> dict[str, str]:
    sample_index = column_index(header, sample_column, "metadata sample")
    sex_index = column_index(header, sex_column, "metadata sex")
    lookup: dict[str, str] = {}

    for line_number, row in rows:
        sample = row[sample_index].strip()
        if not sample:
            raise EnrichmentError(
                f"metadata CSV row {line_number} has an empty sample ID"
            )
        if sample in lookup:
            raise EnrichmentError(
                f"duplicate metadata sample ID '{sample}' at row {line_number}"
            )
        lookup[sample] = normalize_sex(
            row[sex_index],
            f"metadata sample '{sample}' at row {line_number}",
        )
    return lookup


def enrich_manifest(
    header: list[str],
    rows: Iterable[tuple[int, list[str]]],
    sex_lookup: dict[str, str],
    sample_column: int,
    suffix: str,
) -> tuple[list[str], list[list[str]]]:
    sample_index = column_index(header, sample_column, "manifest sample")
    sex_columns = [
        index for index, name in enumerate(header) if name.strip().lower() == "sex"
    ]
    if len(sex_columns) > 1:
        raise EnrichmentError("manifest CSV contains more than one 'sex' column")

    output_header = list(header)
    if sex_columns:
        sex_index = sex_columns[0]
    else:
        sex_index = len(output_header)
        output_header.append("sex")

    output_rows: list[list[str]] = []
    normalized_ids: dict[str, int] = {}
    for line_number, original_row in rows:
        row = list(original_row)
        raw_sample = row[sample_index]
        match_id = normalize_manifest_sample(raw_sample, suffix)
        if not match_id:
            raise EnrichmentError(
                f"manifest CSV row {line_number} has an empty sample ID"
            )
        if match_id in normalized_ids:
            raise EnrichmentError(
                f"manifest rows {normalized_ids[match_id]} and {line_number} "
                f"both normalize to metadata ID '{match_id}'"
            )
        normalized_ids[match_id] = line_number

        if match_id not in sex_lookup:
            raise EnrichmentError(
                f"no metadata match for manifest sample '{raw_sample.strip()}' "
                f"(normalized ID '{match_id}') at row {line_number}"
            )
        metadata_sex = sex_lookup[match_id]

        if sex_index == len(row):
            row.append(metadata_sex)
        else:
            existing_sex = row[sex_index].strip()
            if existing_sex:
                normalized_existing = normalize_sex(
                    existing_sex,
                    f"manifest sample '{raw_sample.strip()}' at row {line_number}",
                )
                if normalized_existing != metadata_sex:
                    raise EnrichmentError(
                        f"conflicting sex for manifest sample "
                        f"'{raw_sample.strip()}': manifest has '{existing_sex}', "
                        f"metadata has '{metadata_sex}'"
                    )
            row[sex_index] = metadata_sex
        output_rows.append(row)

    if not output_rows:
        raise EnrichmentError("manifest CSV has no sample rows")
    return output_header, output_rows


def validate_paths(manifest: Path, metadata: Path, output: Path) -> None:
    output_resolved = output.resolve()
    if output_resolved in {manifest.resolve(), metadata.resolve()}:
        raise EnrichmentError("output path must differ from both input CSV paths")
    if output.exists():
        raise EnrichmentError(f"output path already exists: {output}")
    if not output.parent.is_dir():
        raise EnrichmentError(f"output directory does not exist: {output.parent}")


def run(args: argparse.Namespace) -> int:
    validate_paths(args.manifest, args.metadata, args.output)
    manifest_header, manifest_rows = read_csv(args.manifest, "manifest")
    metadata_header, metadata_rows = read_csv(args.metadata, "metadata")
    sex_lookup = build_sex_lookup(
        metadata_header,
        metadata_rows,
        args.metadata_sample_column,
        args.metadata_sex_column,
    )
    output_header, output_rows = enrich_manifest(
        manifest_header,
        manifest_rows,
        sex_lookup,
        args.manifest_sample_column,
        args.strip_suffix,
    )

    try:
        with args.output.open("x", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, lineterminator="\n")
            writer.writerow(output_header)
            writer.writerows(output_rows)
    except OSError as exc:
        raise EnrichmentError(f"could not write output CSV '{args.output}': {exc}") from exc

    print(f"Wrote {len(output_rows)} enriched sample row(s) to {args.output}")
    return 0


def main(argv: Sequence[str] | None = None) -> int:
    try:
        return run(parse_args(argv))
    except EnrichmentError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
