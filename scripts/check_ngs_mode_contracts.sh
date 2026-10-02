#!/usr/bin/env bash

set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
ngs_workflow="$repo_root/subworkflows/ngs/main.nf"

fail() {
    printf 'ERROR: %s\n' "$*" >&2
    exit 1
}

require_text() {
    local file="$1"
    local text="$2"
    rg --no-ignore -F -q -- "$text" "$file" ||
        fail "missing NGS mode contract in ${file#"$repo_root/"}: $text"
}

reject_text() {
    local file="$1"
    local text="$2"
    if rg --no-ignore -F -q -- "$text" "$file"; then
        fail "obsolete NGS route remains in ${file#"$repo_root/"}: $text"
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
        fail "$label is outside its guarded analysis-mode branch"
    fi
}

[[ -f "$ngs_workflow" ]] || fail 'unified subworkflows/ngs/main.nf is missing'
for obsolete_route in ngs_snp ngs_sv ngs_combined; do
    [[ ! -e "$repo_root/subworkflows/$obsolete_route" ]] ||
        fail "obsolete subworkflow directory remains: subworkflows/$obsolete_route"
done

require_text "$repo_root/main.nf" "include { ALIGNMENT_NGS } from './subworkflows/ngs/main'"
require_text "$repo_root/main.nf" "ngs_analysis_mode = route_settings.effective_mode ?: 'combined'"
reject_text "$repo_root/main.nf" 'ALIGNMENT_ALL_NGS'
reject_text "$repo_root/main.nf" 'ALIGNMENT_NGS_SMALL_VARIANT'
reject_text "$repo_root/main.nf" 'ALIGNMENT_NGS_SV'

require_text "$ngs_workflow" "if (analysis_mode in ['snp', 'combined'])"
require_text "$ngs_workflow" "if (analysis_mode in ['sv', 'combined'])"
require_text "$ngs_workflow" "if (analysis_mode == 'snp')"
require_text "$ngs_workflow" "if (analysis_mode == 'sv')"
require_text "$ngs_workflow" "if (analysis_mode == 'combined')"
require_text "$ngs_workflow" 'small_variant_prio('
require_text "$ngs_workflow" 'sv_prio('
require_text "$ngs_workflow" 'ngs_prio('
require_text "$ngs_workflow" 'variant_html_report_with_repeat('
require_text "$ngs_workflow" 'variant_html_report_with_repeat_and_mito('
require_text "$ngs_workflow" 'failOnMismatch: true'
require_text "$ngs_workflow" 'failOnDuplicate: true'

# Verify that mode-specific calls remain lexically inside their guarded blocks.
snp_start="$(line_number "$ngs_workflow" "if (analysis_mode in ['snp', 'combined'])")"
sv_start="$(line_number "$ngs_workflow" "if (analysis_mode in ['sv', 'combined'])")"
combined_start="$(line_number "$ngs_workflow" "if (analysis_mode == 'combined')")"
require_between 'DeepVariant' "$(line_number "$ngs_workflow" '            DEEPVARIANT_RUNDEEPVARIANT(')" "$snp_start" "$sv_start"
require_between 'SNP prioritization' "$(line_number "$ngs_workflow" '            small_variant_prio(')" "$snp_start" "$sv_start"
require_between 'Manta' "$(line_number "$ngs_workflow" '        MANTA_GERMLINE(')" "$sv_start" "$combined_start"
require_between 'SV prioritization' "$(line_number "$ngs_workflow" '            sv_prio(')" "$sv_start" "$combined_start"
require_between 'combined prioritization' "$(line_number "$ngs_workflow" '        ngs_prio(')" "$combined_start" 1000000
require_between 'combined report' "$(line_number "$ngs_workflow" '            variant_html_report_with_repeat_and_mito(')" "$combined_start" 1000000

printf 'NGS mode contracts OK: one dispatcher preserves SNP, SV, combined, and report branches.\n'
