#!/usr/bin/env bash

set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
status=0

fail() {
    printf 'ERROR: %s\n' "$*" >&2
    status=1
}

count_definitions() {
    local pattern="$1"
    local search_root="$2"
    {
        rg --no-ignore --glob '*.nf' --only-matching "$pattern" "$search_root" || true
    } | wc -l | tr -d '[:space:]'
}

# Check the current tree rather than pinning a historical directory inventory.
process_count="$(count_definitions '^[[:space:]]*process[[:space:]]+[A-Za-z_][A-Za-z0-9_]*[[:space:]]*\{' "$repo_root/modules")"
workflow_count="$(count_definitions '^[[:space:]]*workflow[[:space:]]+[A-Za-z_][A-Za-z0-9_]*[[:space:]]*\{' "$repo_root/subworkflows")"

# Every active DSL2 include must resolve to an existing source file. Nextflow
# accepts a direct file, an extensionless .nf path, or a directory main.nf.
while IFS= read -r include_record; do
    source_file="${include_record%%:*}"
    include_line="${include_record#*:}"
    include_line="${include_line#*:}"
    include_target="$(sed -nE "s/.*from[[:space:]]+['\"]([^'\"]+)['\"].*/\1/p" <<< "$include_line")"
    [[ -n "$include_target" ]] || continue
    resolved_base="$(dirname "$source_file")/$include_target"
    if [[ ! -f "$resolved_base" && ! -f "${resolved_base}.nf" && ! -f "$resolved_base/main.nf" ]]; then
        fail "include target does not resolve: ${source_file#"$repo_root/"} -> $include_target"
    fi
done < <(
    rg --no-ignore -n --glob '*.nf' \
        "^[[:space:]]*include[[:space:]]+.*from[[:space:]]+['\"][^'\"]+['\"]" \
        "$repo_root/main.nf" "$repo_root/modules" "$repo_root/subworkflows"
)

if (( status != 0 )); then
    exit "$status"
fi

printf 'Structure OK: all includes resolve; %s processes / %s subworkflows.\n' "$process_count" "$workflow_count"
