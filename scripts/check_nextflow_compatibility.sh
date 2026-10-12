#!/usr/bin/env bash

set -uo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
nextflow_bin="${NEXTFLOW_BIN:-nextflow}"
compat_root="${COMPAT_ROOT:-/tmp/pipevar-nextflow-compat.$(date +%Y%m%d-%H%M%S)}"
lanes_text="${NEXTFLOW_COMPATIBILITY_LANES:-24.10.6@default 25.04.8@default 25.10.7@v1 25.10.7@v2 26.04.6@default 26.04.6@v1}"
read -r -a lanes <<< "$lanes_text"

if command -v "$nextflow_bin" >/dev/null 2>&1; then
    nextflow_command=("$(command -v "$nextflow_bin")")
elif [[ -x "$nextflow_bin" ]]; then
    nextflow_command=("$nextflow_bin")
elif [[ -f "$nextflow_bin" ]]; then
    nextflow_command=(bash "$nextflow_bin")
else
    printf 'ERROR: Nextflow launcher not found: %s\n' "$nextflow_bin" >&2
    exit 2
fi

java_command="${JAVA_HOME:+$JAVA_HOME/bin/}java"
if ! command -v "$java_command" >/dev/null 2>&1 && [[ ! -x "$java_command" ]]; then
    printf 'ERROR: Java executable not found: %s\n' "$java_command" >&2
    exit 2
fi

java_spec="$($java_command -XshowSettings:properties -version 2>&1 | awk -F'= ' '/java.specification.version =/ {print $2; exit}')"
java_major="${java_spec#1.}"
java_major="${java_major%%.*}"
if [[ ! "$java_major" =~ ^[0-9]+$ ]] || (( java_major < 17 )); then
    printf 'ERROR: PipeVar compatibility checks require Java 17 or newer; detected %s.\n' "${java_spec:-unknown}" >&2
    exit 2
fi

mkdir -p "$compat_root"
summary="$compat_root/summary.tsv"
printf 'lane\tversion\tparser\tcheck\tstatus\texit_code\tlog\n' > "$summary"
overall_status=0

run_logged() {
    local lane="$1"
    local version="$2"
    local parser="$3"
    local check_name="$4"
    local log_file="$5"
    shift 5
    local exit_code

    {
        printf 'COMMAND:'
        printf ' %q' "$@"
        printf '\n'
    } > "$log_file"

    if "$@" >> "$log_file" 2>&1; then
        exit_code=0
        printf '%s\t%s\t%s\t%s\tPASS\t0\t%s\n' "$lane" "$version" "$parser" "$check_name" "$log_file" >> "$summary"
    else
        exit_code=$?
        overall_status=1
        printf '%s\t%s\t%s\t%s\tFAIL\t%s\t%s\n' "$lane" "$version" "$parser" "$check_name" "$exit_code" "$log_file" >> "$summary"
    fi
}

run_logged_expect_failure() {
    local lane="$1"
    local version="$2"
    local parser="$3"
    local check_name="$4"
    local log_file="$5"
    shift 5
    local exit_code

    {
        printf 'COMMAND:'
        printf ' %q' "$@"
        printf '\n'
    } > "$log_file"

    if "$@" >> "$log_file" 2>&1; then
        overall_status=1
        printf '%s\t%s\t%s\t%s\tFAIL\t0\t%s\n' "$lane" "$version" "$parser" "$check_name" "$log_file" >> "$summary"
    else
        exit_code=$?
        printf '%s\t%s\t%s\t%s\tPASS\t%s\t%s\n' "$lane" "$version" "$parser" "$check_name" "$exit_code" "$log_file" >> "$summary"
    fi
}

record_assertion() {
    local lane="$1"
    local version="$2"
    local parser="$3"
    local check_name="$4"
    local log_file="$5"
    local status="$6"

    if [[ "$status" == 'PASS' ]]; then
        printf '%s\t%s\t%s\t%s\tPASS\t0\t%s\n' "$lane" "$version" "$parser" "$check_name" "$log_file" >> "$summary"
    else
        overall_status=1
        printf '%s\t%s\t%s\t%s\tFAIL\t1\t%s\n' "$lane" "$version" "$parser" "$check_name" "$log_file" >> "$summary"
    fi
}

pipevar_revision="$(git -c safe.directory="$repo_root" -C "$repo_root" rev-parse HEAD 2>/dev/null || printf 'unavailable')"
git -c safe.directory="$repo_root" -C "$repo_root" status --short > "$compat_root/source-status.txt" 2>/dev/null || true
printf 'pipevar_revision=%s\n' "$pipevar_revision" > "$compat_root/source-revision.txt"

for lane_spec in "${lanes[@]}"; do
    if [[ "$lane_spec" != *@* ]]; then
        printf 'ERROR: Invalid compatibility lane %q; expected VERSION@default, VERSION@v1, or VERSION@v2.\n' "$lane_spec" >&2
        exit 2
    fi
    version="${lane_spec%%@*}"
    parser="${lane_spec#*@}"
    if [[ "$parser" != 'default' && "$parser" != 'v1' && "$parser" != 'v2' ]]; then
        printf 'ERROR: Invalid parser %q in lane %q.\n' "$parser" "$lane_spec" >&2
        exit 2
    fi
    lane="${version}-${parser}"
    version_root="$compat_root/$lane"
    nxf_home="$version_root/nxf-home"
    run_dir="$version_root/run"
    logs_dir="$version_root/logs"
    work_dir="$version_root/work"
    mkdir -p "$nxf_home" "$run_dir" "$logs_dir" "$work_dir"

    export NXF_VER="$version"
    export NXF_HOME="$nxf_home"
    if [[ "$parser" == 'default' ]]; then
        unset NXF_SYNTAX_PARSER
        parser_environment='unset'
    else
        export NXF_SYNTAX_PARSER="$parser"
        parser_environment="$parser"
    fi

    (
        printf 'requested_nextflow=%s\n' "$version"
        printf 'nextflow_launcher=%s\n' "${nextflow_command[*]}"
        printf 'java_home=%s\n' "${JAVA_HOME:-system}"
        printf 'syntax_parser=%s\n' "$parser_environment"
        printf 'pipevar_revision=%s\n' "$pipevar_revision"
        "$java_command" -version
    ) > "$version_root/environment.txt" 2>&1

    run_logged "$lane" "$version" "$parser" version "$logs_dir/version.log" \
        bash -c 'cd "$1" && shift && "$@"' _ "$run_dir" \
        "${nextflow_command[@]}" -version

    for profile in standard slurm_singularity local_singularity local_docker; do
        run_logged "$lane" "$version" "$parser" "config_${profile}" "$logs_dir/config_${profile}.log" \
            bash -c 'cd "$1" && shift && "$@"' _ "$run_dir" \
            "${nextflow_command[@]}" -C "$repo_root/nextflow.config" config -profile "$profile"
    done

    strict_lane=false
    if [[ "$parser" == 'v2' || ( "$version" == 26.* && "$parser" == 'default' ) ]]; then
        strict_lane=true
    fi
    if [[ "$strict_lane" == true ]]; then
        run_logged "$lane" "$version" "$parser" lint "$logs_dir/lint.log" \
            bash -c 'cd "$1" && shift && "$@"' _ "$run_dir" \
            "${nextflow_command[@]}" lint "$repo_root/main.nf"
    fi

    run_logged "$lane" "$version" "$parser" graph_preview "$logs_dir/graph_preview.log" \
        bash -c 'cd "$1" && shift && "$@"' _ "$run_dir" \
        "${nextflow_command[@]}" -C "$repo_root/nextflow.config" run "$repo_root/main.nf" \
        -profile local_docker -preview -ansi-log false -work-dir "$work_dir/preview" --help

    run_logged_expect_failure "$lane" "$version" "$parser" help_default "$logs_dir/help_default.log" \
        bash -c 'cd "$1" && shift && "$@"' _ "$run_dir" \
        "${nextflow_command[@]}" -C "$repo_root/nextflow.config" run "$repo_root/main.nf" \
        -profile local_docker -preview -ansi-log false -work-dir "$work_dir/help-default"
    if grep -Fq 'ERROR: No input data specified' "$logs_dir/help_default.log" && ! grep -Fxq '  PipeVar_simplified' "$logs_dir/help_default.log"; then
        record_assertion "$lane" "$version" "$parser" help_default_output "$logs_dir/help_default.log" PASS
    else
        record_assertion "$lane" "$version" "$parser" help_default_output "$logs_dir/help_default.log" FAIL
    fi

    run_logged_expect_failure "$lane" "$version" "$parser" help_false "$logs_dir/help_false.log" \
        bash -c 'cd "$1" && shift && "$@"' _ "$run_dir" \
        "${nextflow_command[@]}" -C "$repo_root/nextflow.config" run "$repo_root/main.nf" \
        -profile local_docker -preview -ansi-log false -work-dir "$work_dir/help-false" --help false
    if grep -Fq 'ERROR: No input data specified' "$logs_dir/help_false.log" && ! grep -Fxq '  PipeVar_simplified' "$logs_dir/help_false.log"; then
        record_assertion "$lane" "$version" "$parser" help_false_output "$logs_dir/help_false.log" PASS
    else
        record_assertion "$lane" "$version" "$parser" help_false_output "$logs_dir/help_false.log" FAIL
    fi

    selector_logs=(
        "$logs_dir/graph_preview.log"
        "$logs_dir/help_default.log"
        "$logs_dir/help_false.log"
    )
    if grep -Fq "There's no process matching config selector" "${selector_logs[@]}"; then
        record_assertion "$lane" "$version" "$parser" process_selector_coverage "$logs_dir/graph_preview.log" FAIL
    else
        record_assertion "$lane" "$version" "$parser" process_selector_coverage "$logs_dir/graph_preview.log" PASS
    fi

    run_logged "$lane" "$version" "$parser" strict_syntax_contracts "$logs_dir/strict_syntax_contracts.log" \
        bash "$repo_root/scripts/check_strict_syntax_contracts.sh"
    run_logged "$lane" "$version" "$parser" structure "$logs_dir/structure.log" \
        bash "$repo_root/scripts/check_structure.sh"
    run_logged "$lane" "$version" "$parser" nfcore_contracts "$logs_dir/nfcore_contracts.log" \
        bash "$repo_root/scripts/check_nfcore_contracts.sh"
    run_logged "$lane" "$version" "$parser" sensitive_contracts "$logs_dir/sensitive_contracts.log" \
        bash "$repo_root/scripts/check_sensitive_contracts.sh"
    run_logged "$lane" "$version" "$parser" ngs_mode_contracts "$logs_dir/ngs_mode_contracts.log" \
        bash "$repo_root/scripts/check_ngs_mode_contracts.sh"
    run_logged "$lane" "$version" "$parser" vcf_mode_contracts "$logs_dir/vcf_mode_contracts.log" \
        bash "$repo_root/scripts/check_vcf_mode_contracts.sh"
    run_logged "$lane" "$version" "$parser" long_mode_contracts "$logs_dir/long_mode_contracts.log" \
        bash "$repo_root/scripts/check_long_mode_contracts.sh"
    run_logged "$lane" "$version" "$parser" preannotated_contracts "$logs_dir/preannotated_contracts.log" \
        bash "$repo_root/scripts/check_preannotated_contracts.sh"
    run_logged "$lane" "$version" "$parser" channel_contracts "$logs_dir/channel_contracts.log" \
        bash "$repo_root/scripts/check_channel_contracts.sh" "$nextflow_bin"
    run_logged "$lane" "$version" "$parser" manifest_contracts "$logs_dir/manifest_contracts.log" \
        bash "$repo_root/scripts/check_manifest_contracts.sh" "$nextflow_bin"

    printf '##fileformat=VCFv4.2\n' > "$run_dir/sample.vcf"
    printf 'HP:0001250\n' > "$run_dir/sample.hpo"
    printf '>chr1\nA\n' > "$run_dir/reference.fa"
    printf 'chr1\t1\t6\t1\t2\n' > "$run_dir/reference.fa.fai"
    smoke_log="$logs_dir/smoke.log"
    run_logged "$lane" "$version" "$parser" smoke "$smoke_log" \
        bash -c 'cd "$1" && shift && "$@"' _ "$run_dir" \
        "${nextflow_command[@]}" -C "$repo_root/nextflow.config" \
        run "$repo_root/scripts/fixtures/nextflow_compatibility_smoke.nf" \
        -profile local_docker -ansi-log false -work-dir "$work_dir/smoke" \
        --vcf "$run_dir/sample.vcf" --ref_fa "$run_dir/reference.fa" \
        --hpo "$run_dir/sample.hpo" --mode snp --out_prefix compatibility_smoke
    if ! grep -Fq 'SMOKE|vcf|false|snp' "$smoke_log"; then
        record_assertion "$lane" "$version" "$parser" smoke_output "$smoke_log" FAIL
    else
        record_assertion "$lane" "$version" "$parser" smoke_output "$smoke_log" PASS
    fi

    find "$nxf_home/framework/$version" -maxdepth 1 -type f -name 'nextflow-*-one.jar' \
        -exec sha256sum {} + > "$version_root/artifact-sha256.txt" 2>/dev/null || true
done

printf 'Compatibility results: %s\n' "$summary"
exit "$overall_status"
