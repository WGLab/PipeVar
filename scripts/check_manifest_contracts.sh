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

fixture_dir="$(mktemp -d /tmp/pipevar-manifest-contracts.XXXXXX)"
trap 'rm -rf -- "$fixture_dir"' EXIT

printf '##fileformat=VCFv4.2\n' > "$fixture_dir/sample.vcf"
printf 'HP:0001250\n' > "$fixture_dir/sample.hpo"
printf '>chr1\nA\n' > "$fixture_dir/reference.fa"
printf 'chr1\t1\t6\t1\t2\n' > "$fixture_dir/reference.fa.fai"
printf 'test SIF fixture\n' > "$fixture_dir/phenotagger.sif"

expect_failure() {
    local label="$1"
    local expected="$2"
    shift 2
    local log_file="$fixture_dir/${label}.log"

    if (
        cd "$fixture_dir"
        "${nextflow_command[@]}" \
            -C "$repo_root/nextflow.config" \
            run "$repo_root/main.nf" \
            -ansi-log false \
            -work-dir "$fixture_dir/work-${label}" \
            "$@"
    ) >"$log_file" 2>&1; then
        printf 'ERROR: %s unexpectedly succeeded\n' "$label" >&2
        return 1
    fi
    if ! grep -Fq -- "$expected" "$log_file"; then
        printf 'ERROR: %s did not report expected text: %s\n' "$label" "$expected" >&2
        sed -n '1,120p' "$log_file" >&2
        return 1
    fi
}

expect_preparation() {
    local label="$1"
    shift
    local log_file="$fixture_dir/${label}.log"

    if ! (
        cd "$fixture_dir"
        "${nextflow_command[@]}" \
            -C "$repo_root/nextflow.config" \
            run "$fixture_dir/input-preparation.nf" \
            -ansi-log false \
            -work-dir "$fixture_dir/work-${label}" \
            "$@"
    ) >"$log_file" 2>&1; then
        printf 'ERROR: %s input preparation failed\n' "$label" >&2
        sed -n '1,160p' "$log_file" >&2
        return 1
    fi

    local expected
    for expected in "${EXPECTED_PREPARATION_OUTPUT[@]}"; do
        if ! grep -Fq -- "$expected" "$log_file"; then
            printf 'ERROR: %s did not emit expected record: %s\n' "$label" "$expected" >&2
            sed -n '1,160p' "$log_file" >&2
            return 1
        fi
    done
}

cat > "$fixture_dir/input-preparation.nf" <<EOF
nextflow.enable.dsl = 2
include { INPUT_PREPARATION } from '$repo_root/subworkflows/input_preparation/main'

workflow {
    prepared = INPUT_PREPARATION()
    prepared.route_settings.view { value ->
        'SETTINGS|' + value.route + '|' + value.batch_input + '|' + value.effective_mode
    }
    prepared.input_vcf.view { value -> 'VCF|' + value[0] }
    prepared.clinical_metadata_by_sample.view { value -> 'META|' + value[0] + '|' + value[1] }
}
EOF

printf 'sample,file_path,note_path\n' > "$fixture_dir/header-only.csv"
expect_failure header_only 'must contain at least one sample row' \
    --input_csv "$fixture_dir/header-only.csv" --bam true

printf 'sample,,note_path\ns1,%s,%s\n' "$fixture_dir/sample.vcf" "$fixture_dir/sample.hpo" > "$fixture_dir/blank-header.csv"
expect_failure blank_header 'blank header name(s) at column(s): 2' \
    --input_csv "$fixture_dir/blank-header.csv" --vcf true --mode snp

printf 'sample,file_path,note_path\ns1,%s,%s\ns1,%s,%s\n' \
    "$fixture_dir/sample.vcf" "$fixture_dir/sample.hpo" \
    "$fixture_dir/sample.vcf" "$fixture_dir/sample.hpo" > "$fixture_dir/duplicate.csv"
expect_failure duplicate_sample "duplicate sample 's1'" \
    --input_csv "$fixture_dir/duplicate.csv" --vcf true --mode snp

printf 'sample,file_path,note_path\nbad/sample,%s,%s\n' "$fixture_dir/sample.vcf" "$fixture_dir/sample.hpo" > "$fixture_dir/unsafe.csv"
expect_failure unsafe_sample 'is not filename-safe' \
    --input_csv "$fixture_dir/unsafe.csv" --vcf true --mode snp

printf 'sample,file_path,note_path\ns1,%s,%s\n' "$fixture_dir/sample.vcf" "$fixture_dir/sample.hpo" > "$fixture_dir/legacy.csv"
expect_failure ambiguous_legacy 'choose exactly one of --bam true or --vcf true' \
    --input_csv "$fixture_dir/legacy.csv" --bam true --vcf true --mode snp

expect_failure missing_vcf_mode 'Raw VCF input requires --mode snp or --mode sv' \
    --vcf "$fixture_dir/sample.vcf" --hpo "$fixture_dir/sample.hpo" --out_prefix s1

expect_failure ambiguous_phenotype 'provide exactly one of --note or --hpo' \
    --vcf "$fixture_dir/sample.vcf" --mode snp --note "$fixture_dir/sample.hpo" \
    --hpo "$fixture_dir/sample.hpo" --out_prefix s1

expect_failure invalid_target "Invalid --target/--targeted 'maybe'" \
    --vcf "$fixture_dir/sample.vcf" --mode snp --hpo "$fixture_dir/sample.hpo" \
    --target maybe --out_prefix s1

# Unified input_kind is authoritative; contradictory mode flags must fail early.
printf 'sample,input_kind,phenotype_path,phenotype_format,vcf_path\ns1,vcf_snv,%s,hpo,%s\n' \
    "$fixture_dir/sample.hpo" "$fixture_dir/sample.vcf" > "$fixture_dir/unified.csv"
expect_failure conflicting_mode '--mode conflicts with unified input_kind' \
    --input_csv "$fixture_dir/unified.csv" --mode sv

expect_failure unified_vcf_mito '--mito yes requires BAM/CRAM input' \
    --input_csv "$fixture_dir/unified.csv" --mito yes

expect_failure single_annotated_mito_without_alignment '--mito yes requires BAM/CRAM input' \
    --annotated_snv yes --annovar_txt "$fixture_dir/sample.vcf" \
    --vcf "$fixture_dir/sample.vcf" --hpo "$fixture_dir/sample.hpo" \
    --mito yes --out_prefix s1

expect_failure targeted_vcf_without_reference 'Reference FASTA missing' \
    --vcf "$fixture_dir/sample.vcf" --hpo "$fixture_dir/sample.hpo" \
    --mode snp --target yes --out_prefix s1

expect_failure fractional_batch_size 'must be a positive integer' \
    --vcf "$fixture_dir/sample.vcf" --note "$fixture_dir/sample.hpo" \
    --mode snp --out_prefix s1 --phenotype_extractor phenogpt2 --GPU yes \
    --phenogpt2_batch_size 1.5

expect_failure fractional_gpu_cpus '--gpu_cpus must be a positive integer' \
    --vcf "$fixture_dir/sample.vcf" --hpo "$fixture_dir/sample.hpo" \
    --mode snp --out_prefix s1 --gpu_cpus 2.5

expect_failure fractional_deepvariant_forks '--deepvariant_max_forks must be a positive integer' \
    --vcf "$fixture_dir/sample.vcf" --hpo "$fixture_dir/sample.hpo" \
    --mode snp --out_prefix s1 --deepvariant_max_forks 1.5

expect_failure invalid_gpu_backend 'Invalid --gpu_backend' \
    --vcf "$fixture_dir/sample.vcf" --hpo "$fixture_dir/sample.hpo" \
    --mode snp --out_prefix s1 --gpu_backend banana

expect_failure relative_phenotagger_sif '--phenotagger_sif must be an absolute path' \
    --vcf "$fixture_dir/sample.vcf" --hpo "$fixture_dir/sample.hpo" \
    --mode snp --out_prefix s1 --phenotagger_sif relative.sif

expect_failure missing_phenogpt2_sif '--phenogpt2_sif must be a pre-existing regular file' \
    --vcf "$fixture_dir/sample.vcf" --hpo "$fixture_dir/sample.hpo" \
    --mode snp --out_prefix s1 --phenogpt2_sif "$fixture_dir/missing.sif"

expect_failure docker_rejects_sif 'may be used only with a Singularity profile' \
    -profile local_docker --help --phenotagger_sif "$fixture_dir/phenotagger.sif"

# Equivalent single, one-row, and multi-row inputs must preserve route and sample keys.
EXPECTED_PREPARATION_OUTPUT=('SETTINGS|vcf|false|snp' 'VCF|single_sample' 'META|single_sample|')
expect_preparation single_vcf \
    --vcf "$fixture_dir/sample.vcf" --ref_fa "$fixture_dir/reference.fa" \
    --hpo "$fixture_dir/sample.hpo" --mode snp --out_prefix single_sample

EXPECTED_PREPARATION_OUTPUT=('SETTINGS|vcf|false|snp' 'VCF|staged_sif' 'META|staged_sif|')
expect_preparation valid_staged_sif \
    --vcf "$fixture_dir/sample.vcf" --ref_fa "$fixture_dir/reference.fa" \
    --hpo "$fixture_dir/sample.hpo" --mode snp --out_prefix staged_sif \
    --phenotagger_sif "$fixture_dir/phenotagger.sif"

# Valid yes/no parameters are canonicalized before downstream modules inspect params.
cat > "$fixture_dir/common-sv-normalization.nf" <<EOF
nextflow.enable.dsl = 2
include { INPUT_PREPARATION } from '$repo_root/subworkflows/input_preparation/main'

workflow {
    prepared = INPUT_PREPARATION()
    prepared.route_settings.view { 'COMMON_SV_FILTER|' + params.common_sv_filter }
}
EOF
normalization_log="$fixture_dir/common-sv-normalization.log"
if ! (
        cd "$fixture_dir"
        "${nextflow_command[@]}" \
            -C "$repo_root/nextflow.config" \
            run "$fixture_dir/common-sv-normalization.nf" \
            -ansi-log false \
            -work-dir "$fixture_dir/work-common-sv-normalization" \
            --vcf "$fixture_dir/sample.vcf" --ref_fa "$fixture_dir/reference.fa" \
            --hpo "$fixture_dir/sample.hpo" --mode snp --out_prefix normalized \
            --common_sv_filter ' YES '
    ) >"$normalization_log" 2>&1; then
    printf 'ERROR: common-SV normalization fixture failed\n' >&2
    sed -n '1,160p' "$normalization_log" >&2
    exit 1
fi
grep -Fq 'COMMON_SV_FILTER|yes' "$normalization_log" || {
    printf 'ERROR: common-SV filter was not normalized before downstream use\n' >&2
    sed -n '1,160p' "$normalization_log" >&2
    exit 1
}

printf 'sample,input_kind,phenotype_path,phenotype_format,vcf_path,age\ncohort_one,vcf_snv,%s,hpo,%s,8\n' \
    "$fixture_dir/sample.hpo" "$fixture_dir/sample.vcf" > "$fixture_dir/unified-one.csv"
EXPECTED_PREPARATION_OUTPUT=('SETTINGS|vcf|true|snp' 'VCF|cohort_one' 'META|cohort_one|8y')
expect_preparation one_row_vcf \
    --input_csv "$fixture_dir/unified-one.csv" --ref_fa "$fixture_dir/reference.fa"

printf 'sample,input_kind,phenotype_path,phenotype_format,vcf_path,age_of_onset\ncohort_a,vcf_snv,%s,hpo,%s,4y\ncohort_b,vcf_snv,%s,hpo,%s,\n' \
    "$fixture_dir/sample.hpo" "$fixture_dir/sample.vcf" \
    "$fixture_dir/sample.hpo" "$fixture_dir/sample.vcf" > "$fixture_dir/unified-multi.csv"
EXPECTED_PREPARATION_OUTPUT=('SETTINGS|vcf|true|snp' 'VCF|cohort_a' 'VCF|cohort_b' 'META|cohort_a|4y' 'META|cohort_b|')
expect_preparation multi_row_vcf \
    --input_csv "$fixture_dir/unified-multi.csv" --ref_fa "$fixture_dir/reference.fa"

# The interactive updater must reject quoted input before touching its output.
printf 'sample,input_kind,phenotype_path\ns1,annotated_snv,"note,with-comma.txt"\n' > "$fixture_dir/quoted.csv"
printf 'existing output\n' > "$fixture_dir/updated.csv"
if printf '3\n%s\n%s\n' "$fixture_dir/quoted.csv" "$fixture_dir/updated.csv" | \
    bash "$repo_root/scripts/generate_input_csv.sh" > "$fixture_dir/updater.log" 2>&1; then
    printf 'ERROR: quoted CSV updater input unexpectedly succeeded\n' >&2
    exit 1
fi
grep -Fq 'quoted CSV input is not supported' "$fixture_dir/updater.log"
[[ "$(cat "$fixture_dir/updated.csv")" == 'existing output' ]] || {
    printf 'ERROR: rejected quoted CSV changed the existing output\n' >&2
    exit 1
}

# Legacy updater input may still carry sex, but the rewritten CSV must not
# perpetuate the retired clinical field.
printf '##fileformat=VCFv4.2\n' > "$fixture_dir/s1.hg38_multianno.vcf"
printf 'sample,input_kind,phenotype_path,phenotype_format,age_of_onset,sex,snv_vcf_path\ns1,annotated_snv,%s,hpo,4y,female,%s\n' \
    "$fixture_dir/sample.hpo" "$fixture_dir/sample.vcf" > "$fixture_dir/legacy-sex.csv"
printf '3\n%s\n%s\n%s\n.hg38_multianno.vcf\nn\n' \
    "$fixture_dir/legacy-sex.csv" "$fixture_dir/legacy-sex-updated.csv" \
    "$fixture_dir" | bash "$repo_root/scripts/generate_input_csv.sh" > "$fixture_dir/legacy-sex.log" 2>&1
grep -Fq 'sample,input_kind,phenotype_path,phenotype_format,age_of_onset,snv_vcf_path,sv_vcf_path' \
    "$fixture_dir/legacy-sex-updated.csv"
if grep -Fq ',sex,' "$fixture_dir/legacy-sex-updated.csv" || grep -Fq ',female,' "$fixture_dir/legacy-sex-updated.csv"; then
    printf 'ERROR: retired sex column/value survived CSV updater\n' >&2
    exit 1
fi

printf 'Manifest boundary contracts OK.\n'
