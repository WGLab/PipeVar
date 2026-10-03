#!/usr/bin/env bash

set -euo pipefail

image='beoungl/docker_test:longphase_0.4.0'

command -v docker >/dev/null 2>&1 || {
    printf 'ERROR: docker is required for the LongPhase runtime smoke test.\n' >&2
    exit 1
}

docker image inspect "$image" >/dev/null
docker run --rm --pull=never \
    --entrypoint /bin/bash \
    "$image" \
    -ueo pipefail -c '
        command -v bash >/dev/null
        command -v bcftools >/dev/null
        command -v awk >/dev/null
        test -x /longphase_linux-x64
        /longphase_linux-x64 phase --help >/dev/null
        /longphase_linux-x64 haplotag --help >/dev/null
    '

printf 'LongPhase runtime contract OK: versioned image, utilities, phase, and haplotag.\n'
