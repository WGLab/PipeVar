#!/usr/bin/env bash

set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"

require_text() {
    local file="$1"
    local text="$2"
    if ! rg --no-ignore -F -q -- "$text" "$repo_root/$file"; then
        printf 'ERROR: missing sensitive contract in %s: %s\n' "$file" "$text" >&2
        exit 1
    fi
}

reject_text() {
    local file="$1"
    local text="$2"
    if rg --no-ignore -F -q -- "$text" "$repo_root/$file"; then
        printf 'ERROR: forbidden sensitive-contract text in %s: %s\n' "$file" "$text" >&2
        exit 1
    fi
}

# One validated GATK bundle is acquired per workflow execution and broadcast to
# the three unmodified nf-core GATK stages for both single and CSV inputs.
require_text modules/local/gatk_resource_bundle/main.nf 'process GATK_RESOURCE_BUNDLE'
require_text modules/local/gatk_resource_bundle/main.nf 'md5sum --check --strict'
require_text modules/local/gatk_resource_bundle/main.nf 'test "\$(stat -c'
require_text modules/local/gatk_resource_bundle/main.nf "= '10950827213'"
require_text subworkflows/gatk_snp_calling/main.nf 'GATK_RESOURCE_BUNDLE()'
require_text subworkflows/gatk_snp_calling/main.nf 'GATK4_HAPLOTYPECALLER('
require_text subworkflows/gatk_snp_calling/main.nf 'GATK4_VARIANTRECALIBRATOR('
require_text subworkflows/gatk_snp_calling/main.nf 'GATK4_APPLYVQSR('
require_text subworkflows/gatk_snp_calling/main.nf '--resource:hapmap,known=false,training=true,truth=true,prior=15.0 hapmap_3.3.hg38.vcf.gz'
require_text subworkflows/gatk_snp_calling/main.nf '--resource:omni,known=false,training=true,truth=false,prior=12.0 1000G_omni2.5.hg38.vcf.gz'
require_text subworkflows/gatk_snp_calling/main.nf '--resource:1000G,known=false,training=true,truth=false,prior=10.0 1000G_phase1.snps.high_confidence.hg38.vcf.gz'
require_text subworkflows/gatk_snp_calling/main.nf '--resource:dbsnp,known=true,training=false,truth=false,prior=2.0 Homo_sapiens_assembly38.dbsnp138.vcf.gz'
require_text nextflow.config "ext.args = '-an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -mode SNP'"
require_text nextflow.config "ext.args = '--ts-filter-level 99.0 -mode SNP'"
require_text subworkflows/ngs/main.nf "include { GATK_SMALL_VARIANT_CALLING } from '../gatk_snp_calling'"
require_text subworkflows/input_preparation/main.nf 'params.common_sv_filter = clean_common_sv_filter'

# Every prioritizer consumes the same age-aware contract. Batch and single
# inputs may differ in how age is sourced, but must not select separate
# prioritization implementations.
for family in ngs_prio snp_prio sv_prio; do
    require_text "modules/$family/main.nf" 'age_of_onset'
    reject_text "modules/$family/main.nf" '--sex'
    reject_text "modules/$family/main.nf" '--de-novo'
    if [[ -e "$repo_root/modules/$family/cohort.nf" ]]; then
        printf 'ERROR: obsolete cohort prioritizer remains: modules/%s/cohort.nf\n' "$family" >&2
        exit 1
    fi
    if ! rg -F -q -- '"${age_of_onset}"' "$repo_root/modules/$family/main.nf" && \
       ! rg -F -q -- '"$age_of_onset"' "$repo_root/modules/$family/main.nf"; then
        printf 'ERROR: %s must pass quoted age_of_onset to its prioritizer\n' "$family" >&2
        exit 1
    fi
done

# No active route may retain family/de novo selectors or sex-dependent
# prioritization after the shared age-only contract is installed.
for forbidden in 'DENOVO_' 'denovo_' '--de-novo' 'proband_keys' 'role_column' 'family_column' 'vcf_sample_column' '--sex'; do
    if rg --no-ignore -F -q --glob '*.nf' --glob '*.config' -- "$forbidden" \
        "$repo_root/main.nf" "$repo_root/nextflow.config" "$repo_root/modules" "$repo_root/subworkflows"; then
        printf 'ERROR: forbidden family/de novo/sex behavior remains in active pipeline: %s\n' "$forbidden" >&2
        exit 1
    fi
done
for forbidden_word in sex pedigree proband; do
    if rg --no-ignore -n -w -q --glob '*.nf' --glob '*.config' "$forbidden_word" \
        "$repo_root/main.nf" "$repo_root/nextflow.config" "$repo_root/modules" "$repo_root/subworkflows"; then
        printf 'ERROR: forbidden family/de novo/sex metadata remains in active pipeline: %s\n' "$forbidden_word" >&2
        exit 1
    fi
done

# LongPhase keeps one age-aware prioritizer and batch-only haplotag publication.
require_text modules/local/longphase_support/main.nf 'process LONGPHASE_PRIORITIZE'
reject_text modules/local/longphase_support/main.nf 'process LONGPHASE_PRIORITIZE_COHORT'
longphase_prioritizer_block="$(sed -n '/process LONGPHASE_PRIORITIZE/,$p' "$repo_root/modules/local/longphase_support/main.nf")"
if ! rg -F -q -- 'age_of_onset' <<< "$longphase_prioritizer_block"; then
    printf 'ERROR: LongPhase prioritization must consume age_of_onset\n' >&2
    exit 1
fi
if rg -F -q -- '--sex' <<< "$longphase_prioritizer_block" || rg -F -q -- '--de-novo' <<< "$longphase_prioritizer_block"; then
    printf 'ERROR: LongPhase prioritization must omit sex and de novo behavior\n' >&2
    exit 1
fi
if ! rg -F -q -- '"${age_of_onset}"' <<< "$longphase_prioritizer_block" && \
   ! rg -F -q -- '"$age_of_onset"' <<< "$longphase_prioritizer_block"; then
    printf 'ERROR: LongPhase must pass quoted age_of_onset to its prioritizer\n' >&2
    exit 1
fi
require_text subworkflows/longphase_processing/main.nf 'LONGPHASE_HAPLOTAG('
require_text nextflow.config 'enabled: { meta.batch_input }'
require_text nextflow.config "type == 'ont' ? '--ont' : type == 'pacbio' ? '--pb' : ''"

# Large phenotype-extraction images may be pre-staged only for Singularity.
require_text nextflow.config 'phenogpt2_sif = null'
require_text nextflow.config 'phenotagger_sif = null'
require_text nextflow.config "withName:'phenogpt2'"
require_text nextflow.config "withName:'phenotagger'"
require_text nextflow.config "engine in ['singularity', 'apptainer']"
require_text nextflow.config "pullTimeout = '2h'"
require_text main.nf 'validateStagedSifEngine(workflow.containerEngine'
require_text subworkflows/input_preparation/main.nf "params.phenogpt2_sif = validateStagedSif"
require_text subworkflows/input_preparation/main.nf "params.phenotagger_sif = validateStagedSif"

# Compatibility-sensitive filenames and audit artifacts are unchanged.
for artifact in root tab vcf; do
    require_text modules/cnvnator/main.nf '${out_prefix}_cnvnator.'"$artifact"
done
require_text modules/common_sv_filter/main.nf '${out_prefix}.common_sv_filtered.vcf'
require_text modules/common_sv_filter/main.nf '${out_prefix}.common_sv_removed.vcf'
require_text modules/common_sv_filter/main.nf '${out_prefix}.common_sv_filter.summary.tsv'
require_text modules/variant_html_report/main.nf '${out_prefix}.variant_html_report.html'

# Optional regions are represented as a staged path or [] in active caller
# routes; the historical string sentinel must not return.
reject_text subworkflows/gatk_snp_calling/main.nf '"null"'
reject_text modules/annovar/main.nf '"null"'
reject_text subworkflows/input_preparation/main.nf '"null"'
reject_text subworkflows/input_preparation/main.nf "'null'"

# The no-BND branch must keep real shell continuations; four source backslashes
# would pass text review but break the generated command.
require_text modules/truvari_shortread_sv_merge/main.nf 'sort_plain_vcf \\'
require_text modules/truvari_shortread_sv_merge/main.nf '${out_prefix}.shortread_sv.truvari_merged.vcf \\'

printf 'Sensitive contracts OK: shared GATK resources, clinical metadata, LongPhase, CNV/SV audits, reports, Truvari, and optional regions.\n'
