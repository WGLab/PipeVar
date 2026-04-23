#!/usr/bin/env python3
import argparse
import csv
import gzip
import shutil
import sys
from pathlib import Path


REQUIRED_TXT_COLUMNS = {
    "Chr",
    "Start",
    "Ref",
    "Alt",
    "Otherinfo4",
    "Otherinfo5",
    "Otherinfo7",
    "Otherinfo8",
    "Otherinfo12",
    "Otherinfo13",
}

ANNOTATION_HINT_COLUMNS = {"Func.refGene", "Gene.refGene", "ExonicFunc.refGene"}


def fail(message: str) -> None:
    print(f"ERROR: {message}", file=sys.stderr)
    raise SystemExit(1)


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open("r")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Validate a pre-annotated ANNOVAR TXT/VCF pair.")
    parser.add_argument("--sample", required=True)
    parser.add_argument("--annovar-txt", required=True, type=Path)
    parser.add_argument("--annovar-vcf", required=True, type=Path)
    parser.add_argument("--validated-txt", required=True, type=Path)
    parser.add_argument("--validated-vcf", required=True, type=Path)
    return parser.parse_args()


def validate_txt(txt_path: Path) -> list[tuple[str, str, str, str]]:
    with open_text(txt_path) as handle:
        reader = csv.reader(handle, delimiter="\t")
        try:
            header = next(reader)
        except StopIteration:
            fail(f"{txt_path} is empty.")

        missing = sorted(REQUIRED_TXT_COLUMNS.difference(header))
        if missing:
            fail(
                f"{txt_path} is missing required ANNOVAR columns: {', '.join(missing)}. "
                "Provide a table_annovar multianno TXT generated with -otherinfo."
            )

        if not ANNOTATION_HINT_COLUMNS.intersection(header):
            fail(
                f"{txt_path} does not look like an ANNOVAR multianno TXT because it is missing "
                "annotation columns such as Func.refGene/Gene.refGene."
            )

        idx = {name: header.index(name) for name in ("Otherinfo4", "Otherinfo5", "Otherinfo7", "Otherinfo8")}
        variant_keys: list[tuple[str, str, str, str]] = []
        for row in reader:
            if not row:
                continue
            if len(row) < len(header):
                fail(f"{txt_path} has a truncated row; expected {len(header)} columns, found {len(row)}.")
            key = (
                row[idx["Otherinfo4"]].strip(),
                row[idx["Otherinfo5"]].strip(),
                row[idx["Otherinfo7"]].strip(),
                row[idx["Otherinfo8"]].strip(),
            )
            if all(key):
                variant_keys.append(key)
            if len(variant_keys) >= 25:
                break

    if not variant_keys:
        fail(f"{txt_path} does not contain any variants with populated Otherinfo4/5/7/8 fields.")
    return variant_keys


def normalize_contig(contig: str) -> str:
    return contig.replace("chr", "", 1).upper() if contig.lower().startswith("chr") else contig.upper()


def normalize_variant_key(contig: str, pos: str, ref: str, alt: str) -> tuple[str, str, str, str]:
    return (normalize_contig(contig), str(pos), ref.upper(), alt.upper())


def validate_vcf(vcf_path: Path, expected_keys: list[tuple[str, str, str, str]]) -> None:
    if vcf_path.suffix == ".gz":
        fail("Annotated-SNV mode currently requires an uncompressed ANNOVAR multianno VCF, not .vcf.gz.")

    header_lines = []
    sample_columns = []
    observed_keys: set[tuple[str, str, str, str]] = set()

    with open_text(vcf_path) as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if line.startswith("##"):
                header_lines.append(line)
                continue
            if line.startswith("#CHROM"):
                fields = line.split("\t")
                sample_columns = fields[9:]
                break
        else:
            fail(f"{vcf_path} is missing the #CHROM header line.")

        if len(sample_columns) != 1:
            fail(
                f"{vcf_path} must contain exactly one sample column for annotated-SNV mode; "
                f"found {len(sample_columns)}."
            )

        annovar_header = any("ANNOVAR" in line.upper() for line in header_lines)
        annovar_info = any(
            token in line for token in ("Func.refGene", "Gene.refGene", "ExonicFunc.refGene", "ANNOVAR_DATE")
        )
        if not (annovar_header or annovar_info):
            fail(
                f"{vcf_path} does not look like an ANNOVAR multianno VCF. "
                "Provide the ANNOVAR-generated multianno VCF paired with the TXT."
            )

        for raw_line in handle:
            if raw_line.startswith("#"):
                continue
            fields = raw_line.rstrip("\n").split("\t")
            if len(fields) < 5:
                continue
            chrom, pos, _vid, ref, alt = fields[:5]
            for alt_allele in alt.split(","):
                observed_keys.add(normalize_variant_key(chrom, pos, ref, alt_allele))

    missing = []
    for key in expected_keys:
        if normalize_variant_key(*key) not in observed_keys:
            missing.append(":".join(key))
        if len(missing) >= 5:
            break
    if missing:
        fail(
            f"{vcf_path} does not match the supplied ANNOVAR TXT. Missing variant(s) from the paired VCF: "
            f"{', '.join(missing)}"
        )


def copy_file(src: Path, dest: Path) -> None:
    shutil.copyfile(src, dest)


def main() -> None:
    args = parse_args()
    if not args.annovar_txt.exists():
        fail(f"ANNOVAR TXT not found: {args.annovar_txt}")
    if not args.annovar_vcf.exists():
        fail(f"ANNOVAR VCF not found: {args.annovar_vcf}")

    if "multianno" not in args.annovar_txt.name or "multianno" not in args.annovar_vcf.name:
        fail("Annotated-SNV mode requires ANNOVAR multianno TXT and VCF filenames containing 'multianno'.")

    expected_keys = validate_txt(args.annovar_txt)
    validate_vcf(args.annovar_vcf, expected_keys)

    copy_file(args.annovar_txt, args.validated_txt)
    copy_file(args.annovar_vcf, args.validated_vcf)


if __name__ == "__main__":
    main()
