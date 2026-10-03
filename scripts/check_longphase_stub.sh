#!/usr/bin/env bash

set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
test_root="$repo_root/tests/longphase_stub"
nextflow_image='nextflow/nextflow:25.04.8'

if command -v nextflow >/dev/null 2>&1; then
    work_root="$(mktemp -d "${TMPDIR:-/tmp}/pipevar-longphase-stub.XXXXXX")"
    trap 'rm -rf -- "$work_root"' EXIT

    for platform in ont pacbio; do
        nextflow -C "$test_root/nextflow.config" \
            run "$test_root/main.nf" \
            -profile longphase_stub \
            -work-dir "$work_root/$platform" \
            -stub-run \
            -ansi-log false \
            --type "$platform"
    done
elif command -v docker >/dev/null 2>&1 && docker image inspect "$nextflow_image" >/dev/null 2>&1; then
    for platform in ont pacbio; do
        docker run --rm \
            --entrypoint /usr/local/bin/nextflow \
            -v "$repo_root:/workspace" \
            -w /workspace \
            "$nextflow_image" \
            -C tests/longphase_stub/nextflow.config \
            run tests/longphase_stub/main.nf \
            -profile longphase_stub \
            -work-dir "/tmp/pipevar-longphase-stub-$platform" \
            -stub-run \
            -ansi-log false \
            --type "$platform"
    done
else
    printf 'ERROR: nextflow or the cached %s image is required for the LongPhase stub contract test.\n' "$nextflow_image" >&2
    exit 1
fi

printf 'LongPhase stub contracts OK: ONT and PacBio keyed outputs.\n'
