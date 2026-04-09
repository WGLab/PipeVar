#!/usr/bin/env python3
import argparse
import csv
import os
from collections import defaultdict
from pathlib import Path

import pysam


PREP_DIR = Path(os.environ.get("MITO_DB_PREPARED", "/opt/mito_db/prepared"))


def load_table(name):
    path = PREP_DIR / name
    mapping = defaultdict(list)
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            key = (row["chrom"], row["pos"], row["ref"], row["alt"])
            mapping[key].append(row)
    return mapping


def norm_chrom(value: str) -> str:
    if value in {"MT", "M", "chrM"}:
        return "chrM"
    return value


def parse_sample_fields(record):
    sample = record.samples[0] if record.samples else {}
    ad = sample.get("AD")
    dp = sample.get("DP", record.info.get("DP", "."))
    af = sample.get("AF", record.info.get("AF", "."))
    if isinstance(af, tuple):
        af = ",".join(str(x) for x in af)
    ref_ad = "."
    alt_ads = []
    if ad:
        ref_ad = str(ad[0])
        alt_ads = [str(x) for x in ad[1:]]
    return {
        "caller_dp": str(dp) if dp is not None else ".",
        "caller_ref_ad": ref_ad,
        "caller_alt_ads": alt_ads,
        "caller_af": str(af) if af is not None else ".",
        "caller_tlod": str(record.info.get("TLOD", ".")),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf", required=True)
    parser.add_argument("--out-prefix", required=True)
    parser.add_argument("--hmtvar-status", default="asset_missing")
    args = parser.parse_args()

    mitimpact = load_table("mitimpact_prepared.tsv")
    tapogee = load_table("tapogee_prepared.tsv")
    mitotip = load_table("mitotip_prepared.tsv")
    mitomap = load_table("mitomap_prepared.tsv")

    in_vcf = pysam.VariantFile(args.vcf)
    plain_vcf_path = f"{args.out_prefix}.mito.annotated.vcf"
    out_vcf = pysam.VariantFile(plain_vcf_path, "w", header=in_vcf.header.copy())
    out_vcf.header.info.add("MT_SUPPORT_SCORE", 1, "Float", "Rule-based mitochondrial support score")
    out_vcf.header.info.add("MT_MITOMAP_STATUS", 1, "String", "MITOMAP support status")
    out_vcf.header.info.add("MT_PRIMARY_GENE", 1, "String", "Primary mitochondrial gene symbol from bundled resources")

    fieldnames = [
        "chrom", "pos", "ref", "alt", "variant_key", "gene_symbol", "feature_type", "variant_class_mito",
        "caller_filter", "caller_qual", "caller_dp", "caller_ref_ad", "caller_alt_ad", "caller_af", "caller_tlod",
        "apogee2_score", "apogee2_probability", "apogee2_pred", "tapogee_score", "tapogee_unbiased_score",
        "mitotip_score", "mitotip_quartile", "mitotip_count", "mitotip_percentage", "mitotip_status",
        "mitomap_entry_count", "mitomap_reported_pathogenic", "mitomap_reported_polymorphism",
        "mitomap_disease_status_summary", "mitomap_disease_name_summary", "mitomap_heteroplasmy_summary",
        "mitomap_homoplasmy_summary", "mitomap_citation_count_max", "mitomap_source_tables", "mitomap_status",
        "embedded_mitomap_status", "embedded_mitomap_clinical_info", "clinvar_clnsig", "hgvs",
        "functional_effect_general", "functional_effect_detailed", "hmtvar_status", "mito_annotation_sources",
    ]

    with open(f"{args.out_prefix}.mito.annotated.tsv", "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()

        for record in in_vcf.fetch():
            if not record.alts:
                continue
            sample_fields = parse_sample_fields(record)
            gene_labels = []
            max_support_score = 0.0
            mitomap_status_value = "."

            for alt_index, alt in enumerate(record.alts):
                key = (norm_chrom(record.contig), str(record.pos), record.ref.upper(), alt.upper())
                mitimpact_rows = mitimpact.get(key, [])
                tapogee_rows = tapogee.get(key, [])
                mitotip_rows = mitotip.get(key, [])
                mitomap_rows = mitomap.get(key, [])

                resource_rows = mitimpact_rows + tapogee_rows
                if not resource_rows:
                    resource_rows = [{"gene_symbol": ".", "feature_type": "other", "source": "."}]

                caller_alt_ad = sample_fields["caller_alt_ads"][alt_index] if alt_index < len(sample_fields["caller_alt_ads"]) else "."
                mt_tip = mitotip_rows[0] if mitotip_rows else {}
                mt_map = mitomap_rows[0] if mitomap_rows else {}

                for resource in resource_rows:
                    gene = resource.get("gene_symbol", ".")
                    if gene != ".":
                        gene_labels.append(gene)
                    variant_class = resource.get("feature_type", "other")
                    sources = sorted({x for x in [resource.get("source", "."), "MITOMAP" if mt_map else ".", "MitoTip" if mt_tip else "."] if x != "."})
                    score = 0.0
                    if (mt_map.get("mitomap_reported_pathogenic", "no") == "yes"):
                        score += 5.0
                    if resource.get("apogee2_pred", ".").lower().startswith("path"):
                        score += 2.0
                    if mt_tip.get("mitotip_quartile", ".") == "Q1":
                        score += 2.0
                    max_support_score = max(max_support_score, score)
                    mitomap_status_value = mt_map.get("mitomap_status", mitomap_status_value)

                    writer.writerow(
                        {
                            "chrom": key[0],
                            "pos": key[1],
                            "ref": key[2],
                            "alt": key[3],
                            "variant_key": "|".join(key),
                            "gene_symbol": gene,
                            "feature_type": resource.get("feature_type", "."),
                            "variant_class_mito": variant_class,
                            "caller_filter": ";".join(record.filter.keys()) if record.filter.keys() else ".",
                            "caller_qual": str(record.qual) if record.qual is not None else ".",
                            "caller_dp": sample_fields["caller_dp"],
                            "caller_ref_ad": sample_fields["caller_ref_ad"],
                            "caller_alt_ad": caller_alt_ad,
                            "caller_af": sample_fields["caller_af"],
                            "caller_tlod": sample_fields["caller_tlod"],
                            "apogee2_score": resource.get("apogee2_score", "."),
                            "apogee2_probability": resource.get("apogee2_probability", "."),
                            "apogee2_pred": resource.get("apogee2_pred", "."),
                            "tapogee_score": resource.get("tapogee_score", "."),
                            "tapogee_unbiased_score": resource.get("tapogee_unbiased_score", "."),
                            "mitotip_score": mt_tip.get("mitotip_score", "."),
                            "mitotip_quartile": mt_tip.get("mitotip_quartile", "."),
                            "mitotip_count": mt_tip.get("mitotip_count", "."),
                            "mitotip_percentage": mt_tip.get("mitotip_percentage", "."),
                            "mitotip_status": mt_tip.get("mitotip_status", "."),
                            "mitomap_entry_count": mt_map.get("mitomap_entry_count", "."),
                            "mitomap_reported_pathogenic": mt_map.get("mitomap_reported_pathogenic", "."),
                            "mitomap_reported_polymorphism": mt_map.get("mitomap_reported_polymorphism", "."),
                            "mitomap_disease_status_summary": mt_map.get("mitomap_disease_status_summary", "."),
                            "mitomap_disease_name_summary": mt_map.get("mitomap_disease_name_summary", "."),
                            "mitomap_heteroplasmy_summary": mt_map.get("mitomap_heteroplasmy_summary", "."),
                            "mitomap_homoplasmy_summary": mt_map.get("mitomap_homoplasmy_summary", "."),
                            "mitomap_citation_count_max": mt_map.get("mitomap_citation_count_max", "."),
                            "mitomap_source_tables": mt_map.get("mitomap_source_tables", "."),
                            "mitomap_status": mt_map.get("mitomap_status", "no_match"),
                            "embedded_mitomap_status": resource.get("embedded_mitomap_status", "."),
                            "embedded_mitomap_clinical_info": resource.get("embedded_mitomap_clinical_info", "."),
                            "clinvar_clnsig": resource.get("clinvar_clnsig", "."),
                            "hgvs": resource.get("hgvs", "."),
                            "functional_effect_general": resource.get("functional_effect_general", "."),
                            "functional_effect_detailed": resource.get("functional_effect_detailed", "."),
                            "hmtvar_status": args.hmtvar_status,
                            "mito_annotation_sources": "|".join(sources) if sources else ".",
                        }
                    )

            record.info["MT_SUPPORT_SCORE"] = round(max_support_score, 3)
            record.info["MT_MITOMAP_STATUS"] = mitomap_status_value
            primary_gene = "|".join(sorted(set(gene_labels))) if gene_labels else "."
            record.info["MT_PRIMARY_GENE"] = primary_gene
            out_vcf.write(record)

    out_vcf.close()
    pysam.tabix_index(plain_vcf_path, preset="vcf", force=True)


if __name__ == "__main__":
    main()
