#!/usr/bin/env bash

set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
vcf_workflow="$repo_root/subworkflows/vcf/main.nf"

fail() {
    printf 'ERROR: %s\n' "$*" >&2
    exit 1
}

require_text() {
    local file="$1"
    local text="$2"
    rg --no-ignore -F -q -- "$text" "$file" ||
        fail "missing VCF mode contract in ${file#"$repo_root/"}: $text"
}

reject_text() {
    local file="$1"
    local text="$2"
    if rg --no-ignore -F -q -- "$text" "$file"; then
        fail "obsolete VCF route remains in ${file#"$repo_root/"}: $text"
    fi
}

line_number() {
    local file="$1"
    local text="$2"
    rg --no-ignore -n -F -- "$text" "$file" | sed -n '1{s/:.*//;p;}'
}

require_between() {
    local label="$1"
    local line="$2"
    local start="$3"
    local end="$4"
    if [[ -z "$line" ]] || (( line <= start || line >= end )); then
        fail "$label is outside its guarded VCF analysis-mode branch"
    fi
}

[[ -f "$vcf_workflow" ]] || fail 'unified subworkflows/vcf/main.nf is missing'
for obsolete_route in vcf_snp vcf_sv; do
    [[ ! -e "$repo_root/subworkflows/$obsolete_route" ]] ||
        fail "obsolete subworkflow directory remains: subworkflows/$obsolete_route"
done

require_text "$repo_root/main.nf" "include { ALIGNMENT_VCF } from './subworkflows/vcf/main'"
reject_text "$repo_root/main.nf" 'ALIGNMENT_VCF_SMALL_VARIANT'
reject_text "$repo_root/main.nf" 'ALIGNMENT_VCF_SV'
require_text "$vcf_workflow" "if (analysis_mode == 'snp')"
require_text "$vcf_workflow" "if (analysis_mode == 'sv')"
require_text "$vcf_workflow" 'small_variant_prio('
require_text "$vcf_workflow" 'sv_prio('
require_text "$vcf_workflow" 'failOnMismatch: true'
require_text "$vcf_workflow" 'failOnDuplicate: true'
reject_text "$vcf_workflow" 'variant_html_report'
reject_text "$vcf_workflow" 'EXPANSIONHUNTER'
reject_text "$vcf_workflow" 'ALIGNMENT_NGS_MITO'

snp_start="$(line_number "$vcf_workflow" "if (analysis_mode == 'snp')")"
sv_start="$(line_number "$vcf_workflow" "if (analysis_mode == 'sv')")"
require_between 'SNP annotation' "$(line_number "$vcf_workflow" '        annovar_for_downstream = annovar(')" "$snp_start" "$sv_start"
require_between 'SNP prioritization' "$(line_number "$vcf_workflow" '        small_variant_prio(')" "$snp_start" "$sv_start"
require_between 'SV annotation' "$(line_number "$vcf_workflow" '        annovar_sv_for_downstream = annovar_sv(')" "$sv_start" 1000000
require_between 'SV prioritization' "$(line_number "$vcf_workflow" '        sv_prio(')" "$sv_start" 1000000

printf 'VCF mode contracts OK: one dispatcher preserves mutually exclusive SNP and SV branches.\n'
