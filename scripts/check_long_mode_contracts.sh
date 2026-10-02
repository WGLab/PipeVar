#!/usr/bin/env bash

set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
long_workflow="$repo_root/subworkflows/long/main.nf"

fail() {
    printf 'ERROR: %s\n' "$*" >&2
    exit 1
}

require_text() {
    local file="$1"
    local text="$2"
    rg --no-ignore -F -q -- "$text" "$file" ||
        fail "missing long-read mode contract in ${file#"$repo_root/"}: $text"
}

reject_text() {
    local file="$1"
    local text="$2"
    if rg --no-ignore -F -q -- "$text" "$file"; then
        fail "obsolete long-read route remains in ${file#"$repo_root/"}: $text"
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
        fail "$label is outside its guarded long-read analysis-mode branch"
    fi
}

[[ -f "$long_workflow" ]] || fail 'unified subworkflows/long/main.nf is missing'
for obsolete_route in long_snp long_sv longphase_combined; do
    [[ ! -e "$repo_root/subworkflows/$obsolete_route" ]] ||
        fail "obsolete subworkflow directory remains: subworkflows/$obsolete_route"
done

require_text "$repo_root/main.nf" "include { ALIGNMENT_LONG } from './subworkflows/long/main'"
require_text "$repo_root/main.nf" "long_analysis_mode = route_settings.effective_mode ?: 'combined'"
reject_text "$repo_root/main.nf" 'ALIGNMENT_LONG_SMALL_VARIANT'
reject_text "$repo_root/main.nf" 'ALIGNMENT_LONG_SV'
reject_text "$repo_root/main.nf" 'ALIGNMENT_ALL_LONGPHASE'
require_text "$long_workflow" "if (analysis_mode in ['snp', 'combined'])"
require_text "$long_workflow" "if (analysis_mode in ['sv', 'combined'])"
require_text "$long_workflow" "if (analysis_mode == 'snp')"
require_text "$long_workflow" "if (analysis_mode == 'sv')"
require_text "$long_workflow" "if (analysis_mode == 'combined')"
require_text "$long_workflow" 'small_variant_prio('
require_text "$long_workflow" 'sv_prio('
require_text "$long_workflow" 'LONGPHASE_PROCESSING('
require_text "$long_workflow" 'variant_html_report_with_repeat('
require_text "$long_workflow" 'variant_html_report_with_repeat_and_mito('
require_text "$long_workflow" 'rankvar_nanocaller('
require_text "$long_workflow" '"${meta.id}.sniffles.vcf"'
require_text "$long_workflow" 'failOnMismatch: true'
require_text "$long_workflow" 'failOnDuplicate: true'

snp_start="$(line_number "$long_workflow" "if (analysis_mode in ['snp', 'combined'])")"
sv_start="$(line_number "$long_workflow" "if (analysis_mode in ['sv', 'combined'])")"
combined_start="$(line_number "$long_workflow" "if (analysis_mode == 'combined')")"
require_between 'Clair3' "$(line_number "$long_workflow" '            small_variant_calls = clair3(')" "$snp_start" "$sv_start"
require_between 'SNP prioritization' "$(line_number "$long_workflow" '            small_variant_prio(')" "$snp_start" "$sv_start"
require_between 'Sniffles' "$(line_number "$long_workflow" '        SNIFFLES(')" "$sv_start" "$combined_start"
require_between 'SV prioritization' "$(line_number "$long_workflow" '            sv_prio(')" "$sv_start" "$combined_start"
require_between 'LongPhase' "$(line_number "$long_workflow" '        longphase_result = LONGPHASE_PROCESSING(')" "$combined_start" 1000000
require_between 'combined report' "$(line_number "$long_workflow" '            variant_html_report_with_repeat_and_mito(')" "$combined_start" 1000000

printf 'Long-read mode contracts OK: one workflow preserves SNP, SV, and combined LongPhase branches.\n'
