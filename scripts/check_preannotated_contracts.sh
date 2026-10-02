#!/usr/bin/env bash

set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
shared="$repo_root/subworkflows/preannotated_shared/main.nf"

fail() {
    printf 'ERROR: %s\n' "$*" >&2
    exit 1
}

require_text() {
    local file="$1"
    local text="$2"
    rg --no-ignore -F -q -- "$text" "$file" ||
        fail "missing pre-annotated contract in ${file#"$repo_root/"}: $text"
}

reject_text() {
    local file="$1"
    local text="$2"
    if rg --no-ignore -F -q -- "$text" "$file"; then
        fail "forbidden pre-annotated contract in ${file#"$repo_root/"}: $text"
    fi
}

[[ -f "$shared" ]] || fail 'shared pre-annotated evidence workflow is missing'
require_text "$shared" 'workflow PREANNOTATED_SNV_EVIDENCE'
require_text "$shared" 'rankscore_preannotated('
require_text "$shared" "workflow PREANNOTATED_IMPORTED_SV_CHAIN"
require_text "$shared" "tuple(sample_id, vcf_file, 'preannotated')"
require_text "$shared" 'workflow PREANNOTATED_COMBINED_EVIDENCE'
require_text "$shared" 'failOnMismatch: true'
require_text "$shared" 'failOnDuplicate: true'

for wrapper in \
    subworkflows/preannotated_snp/core.nf \
    subworkflows/preannotated_sv/main.nf \
    subworkflows/preannotated_combined/main.nf \
    subworkflows/preannotated_combined/called_sv.nf; do
    require_text "$repo_root/$wrapper" 'PREANNOTATED_SNV_EVIDENCE('
    reject_text "$repo_root/$wrapper" ' = null'
done

require_text "$repo_root/subworkflows/preannotated_snp/core.nf" 'small_variant_prio('
reject_text "$repo_root/subworkflows/preannotated_snp/core.nf" 'ngs_prio('
for wrapper in subworkflows/preannotated_sv/main.nf subworkflows/preannotated_combined/main.nf; do
    require_text "$repo_root/$wrapper" 'PREANNOTATED_IMPORTED_SV_CHAIN('
done
for wrapper in subworkflows/preannotated_sv/main.nf subworkflows/preannotated_combined/main.nf subworkflows/preannotated_combined/called_sv.nf; do
    require_text "$repo_root/$wrapper" 'PREANNOTATED_COMBINED_EVIDENCE('
done
require_text "$repo_root/subworkflows/preannotated_sv/main.nf" 'variant_html_report_no_repeat('
require_text "$repo_root/subworkflows/preannotated_combined/main.nf" 'EXPANSIONHUNTER('
reject_text "$repo_root/subworkflows/preannotated_combined/main.nf" 'MANTA_GERMLINE('
require_text "$repo_root/subworkflows/preannotated_combined/called_sv.nf" 'MANTA_GERMLINE('
require_text "$repo_root/subworkflows/preannotated_combined/called_sv.nf" 'Channel.empty()'
require_text "$repo_root/subworkflows/preannotated_combined/called_sv.nf" 'tuple(sample_id, vcf_file, "called")'
reject_text "$repo_root/subworkflows/preannotated_combined/called_sv.nf" 'xtea_vcf = null'

printf 'Pre-annotated contracts OK: shared evidence helpers preserve four distinct route contracts without null paths.\n'
