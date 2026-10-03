#!/usr/bin/env bash

set -euo pipefail

profile=""
image_dir=""
config_file="nextflow.config"

phenogpt2_image="beoungl/docker_test:phenogpt2_0.4"
phenotagger_image="beoungl/docker_test:phenotagger"

usage() {
    cat <<'USAGE'
Usage:
  scripts/prepare_phenotype_images.sh \
    --profile=<standard|slurm_singularity|local_singularity|local_docker> \
    [--image-dir=<path>] [--config=<path>]

Docker profiles pull both Docker tags. Singularity profiles pull validated,
checksum-named SIFs into --image-dir and update phenogpt2_sif and
phenotagger_sif in nextflow.config.
USAGE
}

for arg in "$@"; do
    case "$arg" in
        --profile=*)
            profile="${arg#*=}"
            ;;
        --image-dir=*)
            image_dir="${arg#*=}"
            ;;
        --config=*)
            config_file="${arg#*=}"
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "Unknown option: $arg" >&2
            usage >&2
            exit 1
            ;;
    esac
done

case "$profile" in
    standard|slurm_singularity|local_singularity|local_docker) ;;
    "")
        echo "Error: --profile is required." >&2
        usage >&2
        exit 1
        ;;
    *)
        echo "Error: invalid profile: $profile" >&2
        exit 1
        ;;
esac

if [[ ! -f "$config_file" ]]; then
    echo "Error: Nextflow config not found: $config_file" >&2
    exit 1
fi
config_file="$(cd "$(dirname "$config_file")" && pwd -P)/$(basename "$config_file")"

escape_groovy_single_quote() {
    printf '%s' "$1" | sed -e 's/\\/\\\\/g' -e "s/'/\\\\'/g"
}

update_sif_paths() {
    local phenogpt2_path="$1"
    local phenotagger_path="$2"
    local config_dir tmp_file phenogpt2_value phenotagger_value

    config_dir="$(dirname "$config_file")"
    tmp_file="$(mktemp "${TMPDIR:-$config_dir}/pipevar-sif-config.XXXXXX")"

    if [[ -n "$phenogpt2_path" ]]; then
        phenogpt2_value="'$(escape_groovy_single_quote "$phenogpt2_path")'"
        phenotagger_value="'$(escape_groovy_single_quote "$phenotagger_path")'"
    else
        phenogpt2_value="null"
        phenotagger_value="null"
    fi

    if ! PHENOGPT2_VALUE="$phenogpt2_value" PHENOTAGGER_VALUE="$phenotagger_value" awk '
        BEGIN { phenogpt2_done=0; phenotagger_done=0 }
        {
            if ($0 ~ /^[[:space:]]*phenogpt2_sif[[:space:]]*=/) {
                match($0, /^[[:space:]]*/)
                print substr($0, RSTART, RLENGTH) "phenogpt2_sif = " ENVIRON["PHENOGPT2_VALUE"]
                phenogpt2_done=1
                next
            }
            if ($0 ~ /^[[:space:]]*phenotagger_sif[[:space:]]*=/) {
                match($0, /^[[:space:]]*/)
                print substr($0, RSTART, RLENGTH) "phenotagger_sif = " ENVIRON["PHENOTAGGER_VALUE"]
                phenotagger_done=1
                next
            }
            print
        }
        END {
            if (!phenogpt2_done || !phenotagger_done) exit 2
        }
    ' "$config_file" > "$tmp_file"; then
        rm -f -- "$tmp_file"
        echo "Error: failed to update staged SIF paths in $config_file" >&2
        exit 1
    fi

    mv -- "$tmp_file" "$config_file"
}

if [[ "$profile" == "local_docker" ]]; then
    if ! command -v docker >/dev/null 2>&1; then
        echo "Error: Docker is required for profile local_docker." >&2
        exit 1
    fi

    echo "Pulling PhenoGPT2 Docker image: $phenogpt2_image"
    docker pull "$phenogpt2_image"
    echo "Pulling PhenoTagger Docker image: $phenotagger_image"
    docker pull "$phenotagger_image"
    update_sif_paths "" ""
    echo "Prepared PhenoGPT2 and PhenoTagger Docker images."
    exit 0
fi

if [[ -z "$image_dir" ]]; then
    echo "Error: --image-dir is required for a Singularity profile." >&2
    exit 1
fi

if command -v singularity >/dev/null 2>&1; then
    container_runtime="singularity"
elif command -v apptainer >/dev/null 2>&1; then
    container_runtime="apptainer"
else
    echo "Error: Singularity or Apptainer is required for profile $profile." >&2
    exit 1
fi
if ! command -v sha256sum >/dev/null 2>&1; then
    echo "Error: sha256sum is required to stage phenotype images." >&2
    exit 1
fi

mkdir -p "$image_dir"
image_dir="$(cd "$image_dir" && pwd -P)"
if [[ ! -w "$image_dir" ]]; then
    echo "Error: phenotype image directory is not writable: $image_dir" >&2
    exit 1
fi

validate_sif() {
    local sif_path="$1"
    [[ -f "$sif_path" && -r "$sif_path" ]] || return 1
    "$container_runtime" inspect "$sif_path" >/dev/null
    "$container_runtime" sif list "$sif_path" >/dev/null
}

staged_sif_path=""
stage_sif() {
    local slug="$1"
    local source_image="$2"
    local current_link pull_dir partial hash final link_in_pull

    current_link="${image_dir}/${slug}.current.sif"
    if [[ -e "$current_link" || -L "$current_link" ]]; then
        if ! validate_sif "$current_link"; then
            echo "Error: staged-image pointer is invalid: $current_link" >&2
            exit 1
        fi
        staged_sif_path="$(readlink -f "$current_link")"
        echo "Reusing staged image: $staged_sif_path"
        return
    fi

    pull_dir="$(mktemp -d "${image_dir}/.${slug}.pull.XXXXXX")"
    partial="${pull_dir}/${slug}.sif"

    echo "Pulling $source_image into $image_dir"
    if ! "$container_runtime" pull "$partial" "docker://${source_image}"; then
        rm -f -- "$partial"
        rmdir -- "$pull_dir"
        echo "Error: failed to pull $source_image" >&2
        exit 1
    fi
    if ! validate_sif "$partial"; then
        rm -f -- "$partial"
        rmdir -- "$pull_dir"
        echo "Error: pulled image failed SIF validation: $source_image" >&2
        exit 1
    fi

    hash="$(sha256sum "$partial" | awk '{print $1}')"
    final="${image_dir}/${slug}_sha256-${hash}.sif"
    if [[ -e "$final" ]]; then
        if ! printf '%s  %s\n' "$hash" "$final" | sha256sum --check --status; then
            rm -f -- "$partial"
            rmdir -- "$pull_dir"
            echo "Error: checksum-named SIF does not match its filename: $final" >&2
            exit 1
        fi
        rm -f -- "$partial"
    else
        chmod 0444 "$partial"
        mv -- "$partial" "$final"
    fi

    link_in_pull="${pull_dir}/${slug}.current.sif"
    ln -s "$(basename "$final")" "$link_in_pull"
    mv -- "$link_in_pull" "$current_link"
    rmdir -- "$pull_dir"

    staged_sif_path="$final"
    echo "Staged and validated image: $staged_sif_path"
}

stage_sif "phenogpt2_0.4" "$phenogpt2_image"
phenogpt2_sif="$staged_sif_path"
stage_sif "phenotagger" "$phenotagger_image"
phenotagger_sif="$staged_sif_path"

update_sif_paths "$phenogpt2_sif" "$phenotagger_sif"

echo "Prepared phenotype-extraction SIFs and updated $config_file"
echo "  phenogpt2_sif  = $phenogpt2_sif"
echo "  phenotagger_sif = $phenotagger_sif"
