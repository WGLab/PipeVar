#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

docker build -t pipevar_mito/mito_mutect2:0.1.0 "${ROOT_DIR}/docker/mito_mutect2"
docker build -t pipevar_mito/mito_annotation:0.1.0 "${ROOT_DIR}/docker/mito_annotation"
