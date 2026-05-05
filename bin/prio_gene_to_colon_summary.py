#!/usr/bin/env python3
import argparse
import csv
import gzip
from pathlib import Path
import re
import sys


EXPECTED_COLUMNS = {
    "#RANK",
    "GENE",
    "PRIORITY_TIER",
    "MAX_SCORE",
    "RANKVAR",
    "RANKSCORE",
    "VCF_LINE",
}

CPX_TYPE_CLASSES = {
    "CPX_TYPE_INS_iDEL",
    "CPX_TYPE_INVdel",
    "CPX_TYPE_INVdup",
    "CPX_TYPE_dDUP",
    "CPX_TYPE_dDUP_iDEL",
    "CPX_TYPE_delINV",
    "CPX_TYPE_delINVdel",
    "CPX_TYPE_delINVdup",
    "CPX_TYPE_dupINV",
    "CPX_TYPE_dupINVdel",
    "CPX_TYPE_dupINVdup",
    "CPX_TYPE_piDUP_FR",
    "CPX_TYPE_piDUP_RF",
}

CAGI_CLASSES = {
    "SNV",
    "INDEL",
    "MSNV",
    "MINDEL",
    "TRE",
    "BND",
    "CNV",
    "CPX",
    "CTX",
    "DEL",
    "DUP",
    "INS",
    "INS_ME",
    "INS_ME_ALU",
    "INS_ME_LINE1",
    "INS_ME_SVA",
    "INS_UNK",
    "INV",
} | CPX_TYPE_CLASSES

MITO_CONTIGS = {"M", "MT", "CHRM", "CHRMT"}


def open_text(path, mode):
    if str(path).endswith(".gz"):
        return gzip.open(path, mode + "t", encoding="utf-8", newline="")
    return open(path, mode, encoding="utf-8", newline="")


def sanitize(value, missing="."):
    if value is None:
        return missing
    value = str(value).strip()
    if not value or value == ".":
        return missing
    return value.replace(":", "_")


def parse_info(info_field):
    info = {}
    if not info_field or info_field == ".":
        return info
    for item in info_field.split(";"):
        if not item:
            continue
        if "=" in item:
            key, value = item.split("=", 1)
            info[key] = value
        else:
            info[item] = True
    return info


def parse_sample(format_field, sample_field):
    if not format_field or format_field == "." or not sample_field or sample_field == ".":
        return {}
    keys = format_field.split(":")
    values = sample_field.split(":")
    if len(values) < len(keys):
        values.extend(["."] * (len(keys) - len(values)))
    return dict(zip(keys, values))


def first_present(mapping, keys):
    for key in keys:
        value = mapping.get(key)
        if value not in (None, "", "."):
            return value
    return None


def derive_vaf_from_gt(gt, missing="."):
    if gt in ("0/1", "1/0", "0|1", "1|0"):
        return "0.5"
    if gt in ("1/1", "1|1", "1"):
        return "1.0"
    return missing


def normalize_token(value):
    return str(value or "").strip().strip("<>").upper()


def symbolic_alt_token(alt):
    if alt.startswith("<") and alt.endswith(">"):
        return normalize_token(alt)
    return ""


def is_mito_contig(chrom):
    return normalize_token(chrom) in MITO_CONTIGS


def format_chrom(chrom):
    chrom = str(chrom or "").strip()
    if not chrom:
        return chrom
    if chrom.lower().startswith("chr"):
        return chrom
    if normalize_token(chrom) in {"M", "MT"}:
        return "chrM"
    return f"chr{chrom}"


def is_breakend_alt(alt):
    return "[" in alt or "]" in alt


def is_repeat_expansion(alt, info):
    alt_token = symbolic_alt_token(alt)
    if alt_token in {"STR", "TRE", "REPEAT", "REPEAT_EXPANSION", "CNV:TR"}:
        return True
    repeat_keys = {"REPID", "RU", "REPCN", "REPEATID", "REPEAT_ID", "REPEAT_UNIT"}
    if any(key in info for key in repeat_keys):
        return True
    source_text = " ".join(str(info.get(key, "")) for key in ("SOURCE", "SVTYPE", "SO", "TYPE"))
    source_text = source_text.upper()
    return any(marker in source_text for marker in ("STR", "TANDEM_REPEAT", "REPEAT_EXPANSION"))


def normalize_cpx_type(value):
    if value in (None, "", "."):
        return None
    token = str(value).split(",")[0].strip()
    if not token:
        return None
    candidate = token if token.startswith("CPX_TYPE_") else f"CPX_TYPE_{token}"
    return candidate if candidate in CPX_TYPE_CLASSES else None


def insertion_class_from_evidence(alt, info):
    alt_token = symbolic_alt_token(alt)
    sequence_keys = {"SVINSSEQ", "LEFT_SVINSSEQ", "RIGHT_SVINSSEQ", "SEQ", "INSERTED_SEQUENCE"}
    if any(info.get(key) not in (None, "", ".") for key in sequence_keys):
        return "INS"
    evidence = " ".join(
        str(value)
        for key, value in info.items()
        if key.upper()
        in {
            "MEINFO",
            "ME",
            "MEI",
            "MOBILE_ELEMENT",
            "ME_TYPE",
            "INS_TYPE",
            "INSME",
            "SVTYPE",
            "ALT",
        }
    )
    evidence = f"{alt_token} {evidence}".upper()

    has_mei = any(marker in evidence for marker in ("INS:ME", "MOBILE", "MEINFO", "MEI", "ALU", "LINE", "L1", "SVA"))
    if "ALU" in evidence:
        return "INS_ME_ALU"
    if "LINE1" in evidence or "LINE" in evidence or re.search(r"(^|[^A-Z0-9])L1([^A-Z0-9]|$)", evidence):
        return "INS_ME_LINE1"
    if "SVA" in evidence:
        return "INS_ME_SVA"
    if has_mei:
        return "INS_ME"
    if alt_token == "INS":
        return "INS_UNK"
    return "INS"


def sv_class_from_type(svtype, alt, info):
    token = normalize_token(svtype) or symbolic_alt_token(alt)
    token = token.split(":", 1)[0]
    if token == "TRA":
        return "BND"
    if token == "DUP":
        return "DUP"
    if token == "INS":
        return insertion_class_from_evidence(alt, info)
    mapping = {
        "DEL": "DEL",
        "INV": "INV",
        "BND": "BND",
        "CNV": "CNV",
        "CPX": "CPX",
        "CTX": "CTX",
    }
    return mapping.get(token)


def validate_cagi_class(variant_class):
    if variant_class not in CAGI_CLASSES:
        raise ValueError(f"Unsupported CAGI variant class inferred: {variant_class}")
    return variant_class


def classify_variant(chrom, ref, alt, info):
    if "," in alt:
        raise ValueError(f"Unsplit multi-allelic ALT is not supported: {alt}")

    if is_repeat_expansion(alt, info):
        return "TRE"

    cpx_class = normalize_cpx_type(first_present(info, ("CPX_TYPE", "CPXTYPE")))
    if cpx_class:
        return cpx_class

    if is_breakend_alt(alt):
        return "CTX" if normalize_token(info.get("SVTYPE")) == "CTX" else "BND"

    sv_class = sv_class_from_type(info.get("SVTYPE"), alt, info)
    if sv_class:
        return validate_cagi_class(sv_class)

    alt_token = symbolic_alt_token(alt)
    if alt_token:
        if alt_token.startswith("DUP"):
            return "DUP"
        if alt_token.startswith("INS:ME") or alt_token == "INS":
            return insertion_class_from_evidence(alt, info)
        sv_class = sv_class_from_type(alt_token, alt, info)
        if sv_class:
            return validate_cagi_class(sv_class)
        raise ValueError(f"Unsupported symbolic ALT for CAGI class mapping: {alt}")

    if is_mito_contig(chrom):
        return "MSNV" if len(ref) == 1 and len(alt) == 1 else "MINDEL"
    if len(ref) == 1 and len(alt) == 1:
        return "SNV"
    return "INDEL"


def derive_epcr(row, missing="."):
    rank_value = row.get("#RANK")
    try:
        rank = int(str(rank_value).strip())
    except (TypeError, ValueError):
        return missing
    if rank < 1:
        return missing
    score = max(0.0, 1.0 - ((rank - 1) * 0.01))
    return f"{score:.2f}".rstrip("0").rstrip(".") if score != 1.0 else "1.0"


def normalize_header(header):
    return [column.strip() for column in header]


def read_report_rows(input_path):
    with open_text(input_path, "r") as handle:
        first_nonempty = None
        lines = []
        for line in handle:
            if first_nonempty is None and line.strip():
                first_nonempty = line
            lines.append(line)

    if first_nonempty is None:
        return []
    if first_nonempty.startswith("##") or first_nonempty.startswith("#CHROM"):
        raise ValueError(
            "Input looks like a true VCF. Expected PipeVar gene TSV report starting with #RANK."
        )
    if not first_nonempty.startswith("#RANK"):
        raise ValueError("Expected PipeVar gene TSV report starting with #RANK.")

    reader = csv.reader(lines, delimiter="\t")
    try:
        raw_header = next(reader)
    except StopIteration:
        return []
    header = normalize_header(raw_header)
    if not EXPECTED_COLUMNS.issubset(set(header)):
        missing = ", ".join(sorted(EXPECTED_COLUMNS - set(header)))
        raise ValueError(f"Missing expected prio_gene report column(s): {missing}")

    rows = []
    for raw_row in reader:
        if not raw_row or not any(cell.strip() for cell in raw_row):
            continue
        if len(raw_row) >= len(header):
            fixed_row = raw_row[: len(header) - 1]
            fixed_row.append("\t".join(raw_row[len(header) - 1 :]))
            raw_row = fixed_row
        else:
            raw_row.extend([""] * (len(header) - len(raw_row)))
        rows.append({key: value.strip() for key, value in zip(header, raw_row)})
    return rows


def parse_vcf_record(line):
    fields = line.split("\t")
    if len(fields) < 8:
        raise ValueError(f"Embedded VCF record has fewer than 8 columns: {line}")
    chrom, pos, _record_id, ref, alt, _qual, _filter, info_field = fields[:8]
    format_field = fields[8] if len(fields) > 8 else "."
    sample_field = fields[9] if len(fields) > 9 else "."
    return {
        "chrom": chrom,
        "pos": pos,
        "ref": ref,
        "alt": alt,
        "info": parse_info(info_field),
        "format": format_field,
        "sample": parse_sample(format_field, sample_field),
    }


def convert_rows(rows, sample_name, missing="."):
    output_lines = []
    for row in rows:
        vcf_line = row.get("VCF_LINE", "")
        if not vcf_line:
            continue
        for record_line in vcf_line.split(" ||| "):
            record_line = record_line.strip()
            if not record_line:
                continue
            record = parse_vcf_record(record_line)
            info = record["info"]
            sample_map = record["sample"]
            gene = row.get("GENE") or info.get("Gene")
            fields = [
                sanitize(sample_name, missing),
                classify_variant(record["chrom"], record["ref"], record["alt"], info),
                sanitize(format_chrom(record["chrom"]), missing),
                sanitize(record["pos"], missing),
                sanitize(record["ref"], missing),
                sanitize(record["alt"], missing),
                sanitize(sample_map.get("GT"), missing),
                derive_vaf_from_gt(sample_map.get("GT"), missing),
                sanitize(gene, missing),
                derive_epcr(row, missing),
            ]
            output_lines.append(":".join(fields))
    return output_lines


def default_output_path(input_path):
    path = Path(input_path)
    name = path.name
    for suffix in (".prio_gene.vcf.gz", ".prio_gene.vcf"):
        if name.endswith(suffix):
            return str(path.with_name(name[: -len(suffix)] + ".prio_gene.summary.txt"))
    if name.endswith(".gz"):
        return str(path.with_name(name[:-3] + ".summary.txt"))
    return str(path.with_name(name + ".summary.txt"))


def default_sample_name(input_path):
    name = Path(input_path).name
    for suffix in (".prio_gene.vcf.gz", ".prio_gene.vcf", ".vcf.gz", ".vcf", ".gz"):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return Path(input_path).stem


def convert_file(input_path, output_path=None, sample_name=None, missing="."):
    if output_path is None:
        output_path = default_output_path(input_path)
    if sample_name is None:
        sample_name = default_sample_name(input_path)
    rows = read_report_rows(input_path)
    lines = convert_rows(rows, sample_name, missing)
    with open_text(output_path, "w") as handle:
        for line in lines:
            handle.write(line + "\n")
    return output_path


def build_parser():
    parser = argparse.ArgumentParser(
        description="Convert PipeVar *.prio_gene.vcf gene reports to colon-delimited variant summaries."
    )
    parser.add_argument("input", help="PipeVar *.prio_gene.vcf TSV report")
    parser.add_argument("--output", help="Output colon-delimited summary file")
    parser.add_argument("--sample", help="Proband/sample identifier for output")
    parser.add_argument("--missing", default=".", help="Missing value placeholder")
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    try:
        output_path = convert_file(args.input, args.output, args.sample, args.missing)
    except ValueError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    print(f"Wrote {output_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
