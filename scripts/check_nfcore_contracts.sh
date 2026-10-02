#!/usr/bin/env bash

set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"

fail() {
    printf 'ERROR: %s\n' "$*" >&2
    exit 1
}

require_text() {
    local file="$1"
    local text="$2"
    rg --no-ignore -F -q -- "$text" "$repo_root/$file" ||
        fail "missing nf-core contract in $file: $text"
}

expected_revision='94b413835288b0e80bd29256d0a95e4fce020a8b'
revision_count="$({ rg --no-ignore -o -- "$expected_revision" "$repo_root/modules.json" || true; } | wc -l | tr -d '[:space:]')"
[[ "$revision_count" == '9' ]] || fail "modules.json must pin all nine modules to $expected_revision"

expected_nfcore_hash='815af35af5dc79359b62cdfdfea3fb2d22f2a7c56247f846ac0b5a2a7fd76c82'
actual_nfcore_hash="$({
    cd "$repo_root/modules/nf-core"
    find . -type f -printf '%P\0' | sort -z | xargs -0 sha256sum | sha256sum
} | awk '{print $1}')"
[[ "$actual_nfcore_hash" == "$expected_nfcore_hash" ]] ||
    fail 'vendored nf-core files differ from the locked upstream snapshot'

require_text nextflow.config "withName: 'DEEPVARIANT_RUNDEEPVARIANT'"
require_text nextflow.config "ext.args = '--model_type=WGS'"
require_text nextflow.config "pattern: '*.deepvariant.vcf.gz'"
require_text nextflow.config "withName: 'EXPANSIONHUNTER'"
require_text nextflow.config "withName: 'MANTA_GERMLINE'"
require_text nextflow.config "withName: 'SNIFFLES'"
require_text nextflow.config "ext.args = '--allow-overwrite --output-rnames'"
require_text nextflow.config "pattern: '*.{json,vcf}'"
require_text nextflow.config "withName: 'GATK4_HAPLOTYPECALLER'"
require_text nextflow.config "withName: 'GATK4_VARIANTRECALIBRATOR'"
require_text nextflow.config "withName: 'GATK4_APPLYVQSR'"
require_text nextflow.config "withName: 'LONGPHASE_PHASE'"
require_text nextflow.config "withName: 'LONGPHASE_HAPLOTAG'"
require_text nextflow.config "quay.io/biocontainers/longphase:1.7.3--hf5e1c6e_0"

require_text modules/local/legacy_artifact_compat/main.nf 'bgzip --decompress --stdout'
require_text subworkflows/ngs/main.nf 'DEEPVARIANT_RUNDEEPVARIANT.out.vcf.map'
require_text subworkflows/ngs/main.nf '"${meta.id}.json"'
require_text subworkflows/ngs/main.nf '"${meta.id}_manta.vcf"'
require_text subworkflows/long/main.nf '"${meta.id}.sniffles.vcf"'
require_text subworkflows/gatk_snp_calling/main.nf "modules/nf-core/gatk4/haplotypecaller"
require_text subworkflows/gatk_snp_calling/main.nf "modules/nf-core/gatk4/variantrecalibrator"
require_text subworkflows/gatk_snp_calling/main.nf "modules/nf-core/gatk4/applyvqsr"
require_text subworkflows/longphase_processing/main.nf "modules/nf-core/longphase/phase"
require_text subworkflows/longphase_processing/main.nf "modules/nf-core/longphase/haplotag"

printf 'nf-core contracts OK: nine-module revision, vendored integrity, integrations, selectors, and public compatibility artifacts.\n'
