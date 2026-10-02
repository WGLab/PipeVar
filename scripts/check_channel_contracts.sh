#!/usr/bin/env bash

set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
nextflow_bin="${1:-nextflow}"

if command -v "$nextflow_bin" >/dev/null 2>&1 || [[ -x "$nextflow_bin" ]]; then
    nextflow_command=("$nextflow_bin")
elif [[ -f "$nextflow_bin" ]]; then
    nextflow_command=(bash "$nextflow_bin")
else
    printf 'SKIP: Nextflow executable not found: %s\n' "$nextflow_bin"
    exit 0
fi

fixture_dir="$(mktemp -d /tmp/pipevar-channel-contracts.XXXXXX)"
trap 'rm -rf -- "$fixture_dir"' EXIT

touch "$fixture_dir/sample_a.vcf" "$fixture_dir/sample_b.vcf"

cat > "$fixture_dir/channel-contracts.nf" <<'EOF'
nextflow.enable.dsl = 2

workflow {
    if (params.case == 'success') {
        left = Channel.of(
            tuple('sample_b', 'left_b'),
            tuple('sample_a', 'left_a')
        )
        right = Channel.of(
            tuple('sample_a', 'right_a'),
            tuple('sample_b', 'right_b')
        )
        left.join(right, failOnMismatch: true, failOnDuplicate: true)
            .view { sample_id, left_value, right_value ->
                "JOIN|${sample_id}|${left_value}|${right_value}"
            }

        wide_evidence = Channel.of(
            tuple('sample_b', 'rankvar_b', 'rankscore_b', 'clinvar_b', 'phenosv_b', 'sv_b.vcf', 'snv_b.vcf', 'b.hpo', '8y'),
            tuple('sample_a', 'rankvar_a', 'rankscore_a', 'clinvar_a', 'phenosv_a', 'sv_a.vcf', 'snv_a.vcf', 'a.hpo', '4y')
        )
        evidence_context = wide_evidence.map {
            sample_id, rankvar, rankscore, clinvar, phenosv, sv_vcf, snv_vcf, hpo, age ->
            tuple(sample_id, [
                rankvar: rankvar,
                rankscore: rankscore,
                clinvar: clinvar,
                phenosv: phenosv,
                sv_vcf: sv_vcf,
                snv_vcf: snv_vcf,
                hpo: hpo,
                age: age
            ])
        }
        evidence_context.map { sample_id, evidence ->
            tuple(sample_id, evidence.rankvar, evidence.rankscore, evidence.clinvar,
                evidence.phenosv, evidence.sv_vcf, evidence.snv_vcf, evidence.hpo,
                evidence.age)
        }.view { row -> "CONTEXT|${row.join('|')}" }

        collected_vcfs = Channel.of(
            file("${params.fixture_dir}/sample_b.vcf"),
            file("${params.fixture_dir}/sample_a.vcf")
        ).collect()
        bindings = Channel.of(
            [out_prefix: 'sample_a', vcf_path: 'sample_a.vcf'],
            [out_prefix: 'sample_b', vcf_path: 'sample_b.vcf']
        )
        bindings.combine(collected_vcfs).map { combined ->
            def binding = combined[0]
            def all_vcfs = combined[1..-1]
            def vcf_file = all_vcfs.find { it.name == binding.vcf_path }
            if (vcf_file == null) {
                throw new IllegalStateException("missing binding for ${binding.out_prefix}")
            }
            tuple(binding.out_prefix, vcf_file.name)
        }.view { sample_id, filename -> "BINDING|${sample_id}|${filename}" }

    }
    else if (params.case == 'duplicate') {
        Channel.of(tuple('sample_a', 'left'))
            .join(
                Channel.of(tuple('sample_a', 'right_1'), tuple('sample_a', 'right_2')),
                failOnMismatch: true,
                failOnDuplicate: true
            )
            .view()
    }
    else if (params.case == 'missing') {
        Channel.of(tuple('sample_a', 'left'))
            .join(
                Channel.of(tuple('sample_b', 'right')),
                failOnMismatch: true,
                failOnDuplicate: true
            )
            .view()
    }
}
EOF

success_log="$fixture_dir/success.log"
if ! (
        cd "$fixture_dir"
        "${nextflow_command[@]}" run "$fixture_dir/channel-contracts.nf" \
            -ansi-log false \
            -work-dir "$fixture_dir/work-success" \
            --case success \
            --fixture_dir "$fixture_dir"
    ) >"$success_log" 2>&1; then
    printf 'ERROR: channel contract success case failed\n' >&2
    sed -n '1,200p' "$success_log" >&2
    exit 1
fi

for expected in \
    'JOIN|sample_a|left_a|right_a' \
    'JOIN|sample_b|left_b|right_b' \
    'CONTEXT|sample_a|rankvar_a|rankscore_a|clinvar_a|phenosv_a|sv_a.vcf|snv_a.vcf|a.hpo|4y' \
    'CONTEXT|sample_b|rankvar_b|rankscore_b|clinvar_b|phenosv_b|sv_b.vcf|snv_b.vcf|b.hpo|8y' \
    'BINDING|sample_a|sample_a.vcf' \
    'BINDING|sample_b|sample_b.vcf'; do
    if ! grep -Fq -- "$expected" "$success_log"; then
        printf 'ERROR: channel contract output missing: %s\n' "$expected" >&2
        sed -n '1,200p' "$success_log" >&2
        exit 1
    fi
done

for failure_case in duplicate missing; do
    if (
        cd "$fixture_dir"
        "${nextflow_command[@]}" run "$fixture_dir/channel-contracts.nf" \
            -ansi-log false \
            -work-dir "$fixture_dir/work-$failure_case" \
            --case "$failure_case" \
            --fixture_dir "$fixture_dir"
    ) >"$fixture_dir/$failure_case.log" 2>&1; then
        printf 'ERROR: strict %s-key contract unexpectedly succeeded\n' "$failure_case" >&2
        exit 1
    fi
done

printf 'Channel contracts OK: keyed joins, age-only context flattening, and binding reconstruction.\n'
