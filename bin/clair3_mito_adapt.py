#!/usr/bin/env python3
import argparse
import gzip
from collections import OrderedDict


def open_text(path, mode):
    if path.endswith(".gz"):
        return gzip.open(path, mode + "t", encoding="utf-8")
    return open(path, mode, encoding="utf-8")


def parse_meta(line):
    line = line.rstrip("\n")
    key, rest = line[2:].split("=", 1)
    return key, rest


def format_meta(key, rest):
    return f"##{key}={rest}\n"


def parse_mapping(field_names, sample_value):
    keys = field_names.split(":")
    values = sample_value.split(":")
    if len(values) < len(keys):
        values.extend(["."] * (len(keys) - len(values)))
    return OrderedDict(zip(keys, values))


def parse_first_float(value):
    if value in {None, "", "."}:
        return None
    token = value.split(",")[0]
    try:
        return float(token)
    except ValueError:
        return None


def parse_first_int(value):
    if value in {None, "", "."}:
        return None
    token = value.split(",")[0]
    try:
        return int(float(token))
    except ValueError:
        return None


def clamp_non_negative(value):
    return value if value >= 0 else 0


def derive_fields(sample_map):
    dp = parse_first_int(sample_map.get("DP"))
    af = parse_first_float(sample_map.get("AF"))

    ad_value = sample_map.get("AD")
    if ad_value not in {None, "", "."}:
        tokens = ad_value.split(",")
        if len(tokens) >= 2:
            ref_ad = parse_first_int(tokens[0])
            alt_ad = parse_first_int(tokens[1])
            if ref_ad is not None and alt_ad is not None:
                if dp is None:
                    dp = ref_ad + alt_ad
                if af is None and dp and dp > 0:
                    af = alt_ad / float(dp)
                return ref_ad, alt_ad, dp, af

    alt_ad = parse_first_int(sample_map.get("AO"))
    if alt_ad is None:
        alt_ad = parse_first_int(sample_map.get("ALT_DEPTH"))
    ref_ad = parse_first_int(sample_map.get("RO"))
    if ref_ad is None:
        ref_ad = parse_first_int(sample_map.get("REF_DEPTH"))

    if dp is None and ref_ad is not None and alt_ad is not None:
        dp = ref_ad + alt_ad

    if alt_ad is None and af is not None and dp is not None:
        alt_ad = int(round(af * dp))

    if ref_ad is None and dp is not None and alt_ad is not None:
        ref_ad = clamp_non_negative(dp - alt_ad)

    if dp is None and ref_ad is not None and alt_ad is not None:
        dp = ref_ad + alt_ad

    if af is None and dp is not None and dp > 0 and alt_ad is not None:
        af = alt_ad / float(dp)

    if ref_ad is None:
        ref_ad = 0
    if alt_ad is None:
        alt_ad = 0

    return ref_ad, alt_ad, dp, af


def ensure_header_lines(meta_lines):
    required = [
        ("FILTER", '<ID=PASS,Description="All filters passed">'),
        ("FORMAT", '<ID=GT,Number=1,Type=String,Description="Genotype">'),
        ("FORMAT", '<ID=DP,Number=1,Type=Integer,Description="Read depth">'),
        ("FORMAT", '<ID=AF,Number=A,Type=Float,Description="Allele fraction">'),
        ("FORMAT", '<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles">'),
    ]

    seen = {"FILTER": set(), "FORMAT": set()}
    for key, rest in meta_lines:
        if key in seen and rest.startswith("<ID="):
            entry_id = rest[4:].split(",", 1)[0]
            seen[key].add(entry_id)

    updated = list(meta_lines)
    for key, rest in required:
        entry_id = rest[4:].split(",", 1)[0]
        if entry_id not in seen[key]:
            updated.append((key, rest))
    return updated


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    meta_lines = []
    header_line = None
    records = []

    with open_text(args.input, "r") as handle:
        for line in handle:
            if line.startswith("##"):
                meta_lines.append(parse_meta(line))
                continue
            if line.startswith("#CHROM"):
                header_line = line.rstrip("\n")
                continue
            if line.strip():
                records.append(line.rstrip("\n"))

    if header_line is None:
        raise SystemExit("Missing #CHROM header in input VCF")

    meta_lines = ensure_header_lines(meta_lines)

    with open(args.output, "w", encoding="utf-8") as out_handle:
        for key, rest in meta_lines:
            out_handle.write(format_meta(key, rest))
        out_handle.write(header_line + "\n")

        for record in records:
            parts = record.split("\t")
            if len(parts) < 10:
                out_handle.write(record + "\n")
                continue

            format_keys = parts[8]
            sample_map = parse_mapping(format_keys, parts[9])
            ref_ad, alt_ad, dp, af = derive_fields(sample_map)

            if dp is None:
                dp = ref_ad + alt_ad
            if af is None:
                af = (alt_ad / float(dp)) if dp else 0.0

            output_keys = ["GT", "DP", "AF", "AD"]
            updated_sample = {
                "GT": sample_map.get("GT", "./."),
                "DP": str(dp),
                "AF": f"{af:.6f}".rstrip("0").rstrip(".") if af is not None else ".",
                "AD": f"{ref_ad},{alt_ad}",
            }

            parts[8] = ":".join(output_keys)
            parts[9] = ":".join(updated_sample[key] for key in output_keys)
            out_handle.write("\t".join(parts) + "\n")


if __name__ == "__main__":
    main()
