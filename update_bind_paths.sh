#!/bin/bash
set -euo pipefail

annovar_dir=""
phenosv_dir=""
config_file="nextflow.config"

usage() {
    cat <<'USAGE'
Usage:
  ./update_bind_paths.sh --annovar-dir=<path> --phenosv-dir=<path>
  ./update_bind_paths.sh <annovar_dir> <phenosv_dir>

Updates PipeVar local bind-path settings used by Singularity/Docker in
nextflow.config without running the full setup.

Options:
  --annovar-dir=<path>      Existing ANNOVAR directory (contains annotate_variation.pl)
  --phenosv-dir=<path>      Existing PhenoSV resource directory
  -h, --help                Show this help message
USAGE
}

resolve_path() {
    local path_in="$1"
    if [[ -d "$path_in" ]]; then
        (cd "$path_in" && pwd -P)
    else
        return 1
    fi
}

escape_groovy_single_quote() {
    printf "%s" "$1" | sed "s/'/\\\\'/g"
}

update_nextflow_bind_paths() {
    local annovar_path="$1"
    local phenosv_path="$2"
    local annovar_escaped phenosv_escaped tmp_file

    if [[ ! -f "$config_file" ]]; then
        echo "Error: nextflow.config not found in current directory: $(pwd)" >&2
        exit 1
    fi

    annovar_escaped="$(escape_groovy_single_quote "$annovar_path")"
    phenosv_escaped="$(escape_groovy_single_quote "$phenosv_path")"

    tmp_file="$(mktemp /tmp/pipevar-bind-config.XXXXXX)"
    if ! awk -v annovar="$annovar_escaped" -v phenosv="$phenosv_escaped" '
        BEGIN { annovar_done=0; phenosv_done=0 }
        {
            if ($0 ~ /^[[:space:]]*annovar_host_path[[:space:]]*=/) {
                match($0, /^[[:space:]]*/)
                indent = substr($0, RSTART, RLENGTH)
                print indent "annovar_host_path = \047" annovar "\047"
                annovar_done=1
                next
            }
            if ($0 ~ /^[[:space:]]*phenosv_host_path[[:space:]]*=/) {
                match($0, /^[[:space:]]*/)
                indent = substr($0, RSTART, RLENGTH)
                print indent "phenosv_host_path = \047" phenosv "\047"
                phenosv_done=1
                next
            }
            print
        }
        END {
            if (!annovar_done || !phenosv_done) exit 2
        }
    ' "$config_file" > "$tmp_file"; then
        rm -f "$tmp_file"
        echo "Error: failed to update annovar_host_path/phenosv_host_path in ${config_file}" >&2
        exit 1
    fi
    mv "$tmp_file" "$config_file"
}

for arg in "$@"; do
    case "$arg" in
        --annovar-dir=*)
            annovar_dir="${arg#*=}"
            ;;
        --phenosv-dir=*)
            phenosv_dir="${arg#*=}"
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
            if [[ -z "$annovar_dir" ]]; then
                annovar_dir="$arg"
            elif [[ -z "$phenosv_dir" ]]; then
                phenosv_dir="$arg"
            else
                echo "Unexpected extra argument: $arg" >&2
                usage >&2
                exit 1
            fi
            ;;
    esac
done

if [[ -z "$annovar_dir" || -z "$phenosv_dir" ]]; then
    echo "Error: both ANNOVAR and PhenoSV directories are required." >&2
    usage >&2
    exit 1
fi

if [[ ! -d "$annovar_dir" ]]; then
    echo "Error: ANNOVAR directory not found: $annovar_dir" >&2
    exit 1
fi
if [[ ! -f "$annovar_dir/annotate_variation.pl" ]]; then
    echo "Error: annotate_variation.pl not found in: $annovar_dir" >&2
    exit 1
fi
if [[ ! -d "$phenosv_dir" ]]; then
    echo "Error: PhenoSV directory not found: $phenosv_dir" >&2
    exit 1
fi

annovar_dir="$(resolve_path "$annovar_dir")"
phenosv_dir="$(resolve_path "$phenosv_dir")"

update_nextflow_bind_paths "$annovar_dir" "$phenosv_dir"

echo "Updated ${config_file}"
echo "  annovar_host_path        = ${annovar_dir}"
echo "  phenosv_host_path        = ${phenosv_dir}"
echo
echo "These paths will be used by Singularity/Docker bind options in PipeVar profiles."
