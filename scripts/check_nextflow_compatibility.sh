#!/usr/bin/env bash

set -uo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
nextflow_bin="${NEXTFLOW_BIN:-nextflow}"
compat_root="${COMPAT_ROOT:-/tmp/pipevar-nextflow-compat.$(date +%Y%m%d-%H%M%S)}"
versions_text="${NEXTFLOW_VERSIONS:-24.10.6 25.04.8 25.10.7 26.04.6}"
syntax_parser="${NEXTFLOW_SYNTAX_PARSER:-v1}"
read -r -a versions <<< "$versions_text"

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
printf 'version\tcheck\tstatus\texit_code\tlog\n' > "$summary"
overall_status=0
export NXF_SYNTAX_PARSER="$syntax_parser"

run_logged() {
    local version="$1"
    local check_name="$2"
    local log_file="$3"
    shift 3
    local exit_code

    {
        printf 'COMMAND:'
        printf ' %q' "$@"
        printf '\n'
    } > "$log_file"

    if "$@" >> "$log_file" 2>&1; then
        exit_code=0
        printf '%s\t%s\tPASS\t0\t%s\n' "$version" "$check_name" "$log_file" >> "$summary"
    else
        exit_code=$?
        overall_status=1
        printf '%s\t%s\tFAIL\t%s\t%s\n' "$version" "$check_name" "$exit_code" "$log_file" >> "$summary"
    fi
}

for version in "${versions[@]}"; do
    version_root="$compat_root/$version"
    nxf_home="$version_root/nxf-home"
    run_dir="$version_root/run"
    logs_dir="$version_root/logs"
    work_dir="$version_root/work"
    mkdir -p "$nxf_home" "$run_dir" "$logs_dir" "$work_dir"

    export NXF_VER="$version"
    export NXF_HOME="$nxf_home"

    (
        printf 'requested_nextflow=%s\n' "$version"
        printf 'nextflow_launcher=%s\n' "${nextflow_command[*]}"
        printf 'java_home=%s\n' "${JAVA_HOME:-system}"
        printf 'syntax_parser=%s\n' "$NXF_SYNTAX_PARSER"
        "$java_command" -version
    ) > "$version_root/environment.txt" 2>&1

    run_logged "$version" version "$logs_dir/version.log" \
        "${nextflow_command[@]}" -version

    for profile in standard slurm_singularity local_singularity local_docker; do
        run_logged "$version" "config_${profile}" "$logs_dir/config_${profile}.log" \
            bash -c 'cd "$1" && shift && "$@"' _ "$run_dir" \
            "${nextflow_command[@]}" -C "$repo_root/nextflow.config" config -profile "$profile"
    done

    run_logged "$version" graph_preview "$logs_dir/graph_preview.log" \
        bash -c 'cd "$1" && shift && "$@"' _ "$run_dir" \
        "${nextflow_command[@]}" -C "$repo_root/nextflow.config" run "$repo_root/main.nf" \
        -profile local_docker -preview -ansi-log false -work-dir "$work_dir/preview" --help

    run_logged "$version" structure "$logs_dir/structure.log" \
        bash "$repo_root/scripts/check_structure.sh"
    run_logged "$version" nfcore_contracts "$logs_dir/nfcore_contracts.log" \
        bash "$repo_root/scripts/check_nfcore_contracts.sh"
    run_logged "$version" sensitive_contracts "$logs_dir/sensitive_contracts.log" \
        bash "$repo_root/scripts/check_sensitive_contracts.sh"
    run_logged "$version" ngs_mode_contracts "$logs_dir/ngs_mode_contracts.log" \
        bash "$repo_root/scripts/check_ngs_mode_contracts.sh"
    run_logged "$version" vcf_mode_contracts "$logs_dir/vcf_mode_contracts.log" \
        bash "$repo_root/scripts/check_vcf_mode_contracts.sh"
    run_logged "$version" long_mode_contracts "$logs_dir/long_mode_contracts.log" \
        bash "$repo_root/scripts/check_long_mode_contracts.sh"
    run_logged "$version" preannotated_contracts "$logs_dir/preannotated_contracts.log" \
        bash "$repo_root/scripts/check_preannotated_contracts.sh"
    run_logged "$version" channel_contracts "$logs_dir/channel_contracts.log" \
        bash "$repo_root/scripts/check_channel_contracts.sh" "$nextflow_bin"
    run_logged "$version" manifest_contracts "$logs_dir/manifest_contracts.log" \
        bash "$repo_root/scripts/check_manifest_contracts.sh" "$nextflow_bin"

    printf '##fileformat=VCFv4.2\n' > "$run_dir/sample.vcf"
    printf 'HP:0001250\n' > "$run_dir/sample.hpo"
    printf '>chr1\nA\n' > "$run_dir/reference.fa"
    printf 'chr1\t1\t6\t1\t2\n' > "$run_dir/reference.fa.fai"
    smoke_log="$logs_dir/smoke.log"
    run_logged "$version" smoke "$smoke_log" \
        bash -c 'cd "$1" && shift && "$@"' _ "$run_dir" \
        "${nextflow_command[@]}" -C "$repo_root/nextflow.config" \
        run "$repo_root/scripts/fixtures/nextflow_compatibility_smoke.nf" \
        -profile local_docker -ansi-log false -work-dir "$work_dir/smoke" \
        --vcf "$run_dir/sample.vcf" --ref_fa "$run_dir/reference.fa" \
        --hpo "$run_dir/sample.hpo" --mode snp --out_prefix compatibility_smoke
    if ! grep -Fq 'SMOKE|vcf|false|snp' "$smoke_log"; then
        overall_status=1
        printf '%s\tsmoke_output\tFAIL\t1\t%s\n' "$version" "$smoke_log" >> "$summary"
    else
        printf '%s\tsmoke_output\tPASS\t0\t%s\n' "$version" "$smoke_log" >> "$summary"
    fi

    find "$nxf_home/framework/$version" -maxdepth 1 -type f -name 'nextflow-*-one.jar' \
        -exec sha256sum {} + > "$version_root/artifact-sha256.txt" 2>/dev/null || true
done

printf 'Compatibility results: %s\n' "$summary"
exit "$overall_status"
