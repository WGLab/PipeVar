#!/usr/bin/env python3
import argparse
import csv
from collections import defaultdict
from pathlib import Path


def norm_chrom(value: str) -> str:
    value = (value or "").strip()
    if value in {"MT", "M", "chrM"}:
        return "chrM"
    return value


def write_tsv(path: Path, fieldnames, rows) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def prep_mitimpact(raw_dir: Path, out_dir: Path) -> None:
    src = raw_dir / "MitImpact_db_3.1.3.txt"
    rows = []
    with src.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            rows.append(
                {
                    "chrom": norm_chrom(row.get("Chr", "")),
                    "pos": row.get("Start", ""),
                    "ref": row.get("Ref", "").upper(),
                    "alt": row.get("Alt", "").upper(),
                    "gene_symbol": row.get("Gene_symbol", "."),
                    "hgvs": row.get("HGVS", "."),
                    "functional_effect_general": row.get("Functional_effect_general", "."),
                    "functional_effect_detailed": row.get("Functional_effect_detailed", "."),
                    "apogee2_score": row.get("APOGEE2_score", "."),
                    "apogee2_probability": row.get("APOGEE2_probability", "."),
                    "apogee2_pred": row.get("APOGEE2", "."),
                    "clinvar_clnsig": row.get("Clinvar_CLNSIG", "."),
                    "embedded_mitomap_status": row.get("MITOMAP_Disease_Status", "."),
                    "embedded_mitomap_clinical_info": row.get("MITOMAP_Disease_Clinical_info", "."),
                    "feature_type": "protein_coding",
                    "source": "MitImpact",
                }
            )
    write_tsv(out_dir / "mitimpact_prepared.tsv", list(rows[0].keys()) if rows else [], rows)


def prep_tapogee(raw_dir: Path, out_dir: Path) -> None:
    src = raw_dir / "t-APOGEE_2024.0.1.txt"
    rows = []
    with src.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            rows.append(
                {
                    "chrom": norm_chrom(row.get("Chr", "")),
                    "pos": row.get("Pos", ""),
                    "ref": row.get("Ref", "").upper(),
                    "alt": row.get("Alt", "").upper(),
                    "gene_symbol": row.get("Gene_symbol", "."),
                    "tapogee_score": row.get("t-APOGEE score", "."),
                    "tapogee_unbiased_score": row.get("t-APOGEE unbiased score", "."),
                    "feature_type": "mt_trna",
                    "source": "t-APOGEE",
                }
            )
    write_tsv(out_dir / "tapogee_prepared.tsv", list(rows[0].keys()) if rows else [], rows)


def prep_mitotip(raw_dir: Path, out_dir: Path) -> None:
    src = raw_dir / "mitotip_scores.txt"
    rows = []
    with src.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            alt = row.get("Alt", "").upper()
            if alt == ":":
                continue
            rows.append(
                {
                    "chrom": "chrM",
                    "pos": row.get("Position", ""),
                    "ref": row.get("rCRS", "").upper(),
                    "alt": alt,
                    "mitotip_score": row.get("MitoTIP_Score", "."),
                    "mitotip_quartile": row.get("Quartile", "."),
                    "mitotip_count": row.get("Count", "."),
                    "mitotip_percentage": row.get("Percentage", "."),
                    "mitotip_status": row.get("Mitomap_Status", "."),
                    "feature_type": "mt_trna",
                    "source": "MitoTip",
                }
            )
    write_tsv(out_dir / "mitotip_prepared.tsv", list(rows[0].keys()) if rows else [], rows)


def prep_mitomap(raw_dir: Path, out_dir: Path) -> None:
    src = raw_dir / "mitomap_variant_evidence.tsv"
    grouped = defaultdict(list)
    with src.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            key = (
                norm_chrom(row.get("chrom", "")),
                row.get("pos", ""),
                row.get("ref", "").upper(),
                row.get("alt", "").upper(),
            )
            grouped[key].append(row)

    rows = []
    for key, values in sorted(grouped.items()):
        disease_names = sorted({v.get("disease_name", "").strip() for v in values if v.get("disease_name", "").strip()})
        disease_status = sorted({v.get("disease_status", "").strip() for v in values if v.get("disease_status", "").strip()})
        heteroplasmy = sorted({v.get("heteroplasmy", "").strip() for v in values if v.get("heteroplasmy", "").strip()})
        homoplasmy = sorted({v.get("homoplasmy", "").strip() for v in values if v.get("homoplasmy", "").strip()})
        source_tables = sorted({v.get("mitomap_source_table", "").strip() for v in values if v.get("mitomap_source_table", "").strip()})
        citation_count_max = max([int(v.get("citation_count", "0") or 0) for v in values] or [0])
        rows.append(
            {
                "chrom": key[0],
                "pos": key[1],
                "ref": key[2],
                "alt": key[3],
                "mitomap_entry_count": str(len(values)),
                "mitomap_reported_pathogenic": "yes" if any((v.get("reported_as_pathogenic", "") or "").lower() == "yes" for v in values) else "no",
                "mitomap_reported_polymorphism": "yes" if any((v.get("reported_as_polymorphism", "") or "").lower() == "yes" for v in values) else "no",
                "mitomap_disease_status_summary": "|".join(disease_status) if disease_status else ".",
                "mitomap_disease_name_summary": " / ".join(disease_names) if disease_names else ".",
                "mitomap_heteroplasmy_summary": "|".join(heteroplasmy) if heteroplasmy else ".",
                "mitomap_homoplasmy_summary": "|".join(homoplasmy) if homoplasmy else ".",
                "mitomap_citation_count_max": str(citation_count_max),
                "mitomap_source_tables": "|".join(source_tables) if source_tables else ".",
                "mitomap_status": "matched" if values else "no_match",
                "source": "MITOMAP",
            }
        )

    write_tsv(out_dir / "mitomap_prepared.tsv", list(rows[0].keys()) if rows else [], rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--raw-dir", required=True)
    parser.add_argument("--out-dir", required=True)
    args = parser.parse_args()

    raw_dir = Path(args.raw_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    prep_mitimpact(raw_dir, out_dir)
    prep_tapogee(raw_dir, out_dir)
    prep_mitotip(raw_dir, out_dir)
    prep_mitomap(raw_dir, out_dir)


if __name__ == "__main__":
    main()
