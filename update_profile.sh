#!/bin/bash
set -euo pipefail

selected_profile=""
config_file="nextflow.config"

usage() {
    cat <<'USAGE'
Usage:
  ./update_profile.sh --profile=<standard|slurm_singularity|local_singularity|local_docker>
  ./update_profile.sh <standard|slurm_singularity|local_singularity|local_docker>

Updates only profiles.standard in nextflow.config so running PipeVar without
-profile uses the selected backend behavior.
USAGE
}

make_tmpfile() {
    local config_dir tmp_base
    config_dir="$(cd "$(dirname "$config_file")" && pwd -P)"
    tmp_base="${TMPDIR:-$config_dir}"
    mkdir -p "$tmp_base"
    mktemp "$tmp_base/pipevar-profile-config.XXXXXX"
}

update_standard_profile() {
    local profile="$1"
    local tmp_file

    if [[ ! -f "$config_file" ]]; then
        echo "Error: nextflow.config not found in current directory: $(pwd)" >&2
        exit 1
    fi

    tmp_file="$(make_tmpfile)"
    if ! awk -v profile="$profile" '
        function emit_standard(indent, inner) {
            inner = indent "    "
            print indent "standard {"
            if (profile == "standard" || profile == "slurm_singularity") {
                print inner "process.executor = \047slurm\047"
                print inner "singularity.enabled = true"
                print inner "singularity.autoMounts = true"
                print inner "singularity.runOptions = \"--bind ${params.annovar_host_path}:/annovar,${params.phenosv_host_path}:/PhenoSV/train_data\""
                print inner "docker.enabled = false"
            } else if (profile == "local_singularity") {
                print inner "process.executor = \047local\047"
                print inner "singularity.enabled = true"
                print inner "singularity.autoMounts = true"
                print inner "singularity.runOptions = \"--bind ${params.annovar_host_path}:/annovar,${params.phenosv_host_path}:/PhenoSV/train_data\""
                print inner "docker.enabled = false"
            } else if (profile == "local_docker") {
                print inner "process.executor = \047local\047"
                print inner "docker.enabled = true"
                print inner "docker.runOptions = \"-v ${params.annovar_host_path}:/annovar -v ${params.phenosv_host_path}:/PhenoSV/train_data\""
                print inner "singularity.enabled = false"
            } else {
                return 1
            }
            print indent "}"
            return 0
        }
        BEGIN { in_standard=0; replaced=0 }
        {
            if (!replaced && $0 ~ /^[[:space:]]*standard[[:space:]]*\{[[:space:]]*$/) {
                match($0, /^[[:space:]]*/)
                indent = substr($0, RSTART, RLENGTH)
                if (emit_standard(indent) != 0) exit 3
                in_standard=1
                replaced=1
                next
            }
            if (in_standard) {
                if ($0 ~ /^[[:space:]]*\}[[:space:]]*$/) {
                    in_standard=0
                }
                next
            }
            print
        }
        END {
            if (!replaced) exit 4
        }
    ' "$config_file" > "$tmp_file"; then
        rm -f "$tmp_file"
        echo "Error: failed to update profiles.standard in ${config_file}" >&2
        exit 1
    fi
    mv "$tmp_file" "$config_file"
}

for arg in "$@"; do
    case "$arg" in
        --profile=*)
            selected_profile="${arg#*=}"
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        --*)
            echo "Unknown option: $arg" >&2
            usage >&2
            exit 1
            ;;
        *)
            if [[ -z "$selected_profile" ]]; then
                selected_profile="$arg"
            else
                echo "Unexpected extra argument: $arg" >&2
                usage >&2
                exit 1
            fi
            ;;
    esac
done

case "$selected_profile" in
    standard|slurm_singularity|local_singularity|local_docker) ;;
    "")
        echo "Error: profile is required." >&2
        usage >&2
        exit 1
        ;;
    *)
        echo "Error: invalid profile: $selected_profile" >&2
        usage >&2
        exit 1
        ;;
esac

update_standard_profile "$selected_profile"

echo "Updated ${config_file}"
echo "  profiles.standard backend = ${selected_profile}"
echo
echo "Running PipeVar without -profile will now use this backend behavior."
