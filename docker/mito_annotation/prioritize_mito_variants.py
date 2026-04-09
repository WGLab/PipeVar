#!/usr/bin/env python3
import argparse
import csv


def parse_float(value, default=0.0):
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def parse_int(value, default=0):
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return default


def score_row(row):
    score = 0.0
    if row.get("mitomap_reported_pathogenic") == "yes":
        score += 5.0
    if (row.get("apogee2_pred", "") or "").lower().startswith("path"):
        score += 2.0
    if row.get("mitotip_quartile") == "Q1":
        score += 2.0
    if row.get("clinvar_clnsig", ".") not in {".", ""}:
        score += 1.0
    score += min(parse_float(row.get("caller_af")) * 2.0, 2.0)
    return round(score, 3)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--min-vaf", type=float, required=True)
    parser.add_argument("--min-depth", type=int, required=True)
    parser.add_argument("--min-alt-reads", type=int, required=True)
    args = parser.parse_args()

    rows = []
    with open(args.input, "r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = reader.fieldnames or []
        for row in reader:
            if parse_float(row.get("caller_af")) < args.min_vaf:
                continue
            if parse_int(row.get("caller_dp")) < args.min_depth:
                continue
            if parse_int(row.get("caller_alt_ad")) < args.min_alt_reads:
                continue
            row["mt_support_score"] = score_row(row)
            if parse_float(row["mt_support_score"]) >= 7:
                row["mt_support_class"] = "high"
            elif parse_float(row["mt_support_score"]) >= 4:
                row["mt_support_class"] = "medium"
            else:
                row["mt_support_class"] = "low"
            rows.append(row)

    rows.sort(
        key=lambda r: (
            -parse_float(r.get("mt_support_score")),
            -parse_float(r.get("caller_af")),
            -parse_int(r.get("caller_dp")),
            r.get("variant_key", ""),
        )
    )

    output_fields = (fieldnames or []) + ["mt_support_score", "mt_support_class"]
    with open(args.out, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=output_fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
