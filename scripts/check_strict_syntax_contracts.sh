#!/usr/bin/env bash

set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
nextflow_sources=(
    "$repo_root/main.nf"
    "$repo_root/modules"
    "$repo_root/subworkflows"
    "$repo_root/scripts/fixtures"
    "$repo_root/tests"
)

fail_matches() {
    local description="$1"
    local pattern="$2"
    local matches

    matches="$(rg --glob '*.nf' --line-number --pcre2 "$pattern" "${nextflow_sources[@]}" || true)"
    if [[ -n "$matches" ]]; then
        printf 'ERROR: %s\n%s\n' "$description" "$matches" >&2
        exit 1
    fi
}

fail_matches 'Deprecated uppercase Channel factory remains.' '(?s)\bChannel\s*\.'
fail_matches 'Implicit closure parameter `it` remains.' '(?s)\{\s*it\b'
fail_matches 'Java/Groovy-style typed closure parameter remains.' '\{[^\n{}]*\b(?:String|Boolean|boolean|Integer|Long|Double|Float|BigDecimal|File|Path|Map|List|Set|Closure)\s+[A-Za-z_]\w*\s*(?:,|->)'
fail_matches 'Undeclared `prefix` assignment remains in Nextflow code.' '^\s*prefix\s*='
fail_matches 'Java/Groovy-style local variable declaration remains.' '^\s*(?:byte|short|int|long|float|double|boolean|char|String)\s+[A-Za-z_]\w*\s*(?:=|$)'

if rg --line-number --pcre2 '\bas\s+(?:byte|short|int|long|float|double|boolean|char|String|Integer|Long|Float|Double|Boolean)\b' "$repo_root/nextflow.config"; then
    printf 'ERROR: Unsupported soft cast remains in nextflow.config.\n' >&2
    exit 1
fi

printf 'Strict-syntax source contracts OK: declarations, closure parameters, and channel factories are dual-parser compatible.\n'
