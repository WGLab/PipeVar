#!/bin/bash
set -euo pipefail

version=""
non_interactive="no"
selected_profile=""
annovar_dir=""
phenosv_dir=""
annovar_bind_path=""
phenosv_bind_path=""
config_file=""

usage() {
    cat <<'USAGE'
Usage:
  ./setup.sh [light] [options]

Options:
  --non-interactive
  --profile=<standard|slurm_singularity|local_singularity|local_docker>
  --annovar-dir=<path>      ANNOVAR install directory (contains annotate_variation.pl)
  --phenosv-dir=<path>      Directory to download/extract PhenoSV resources
  --annovar-bind=<path>     Optional override bind source for /annovar
  --phenosv-bind=<path>     Optional override bind source for /PhenoSV/train_data
USAGE
}

escape_groovy_single_quote() {
    printf "%s" "$1" | sed "s/'/\\\\'/g"
}

make_tmpfile() {
    local config_dir tmp_base
    config_dir="$(cd "$(dirname "$config_file")" && pwd -P)"
    tmp_base="${TMPDIR:-$config_dir}"
    mkdir -p "$tmp_base"
    mktemp "$tmp_base/pipevar-nextflow-config.XXXXXX"
}

update_nextflow_config() {
    local profile="$1"
    local annovar_path="$2"
    local phenosv_path="$3"
    local annovar_escaped phenosv_escaped tmp_file

    config_file="${current_folder}/nextflow.config"
    if [[ ! -f "$config_file" ]]; then
        echo "Error: nextflow.config not found: $config_file" >&2
        exit 1
    fi

    annovar_escaped="$(escape_groovy_single_quote "$annovar_path")"
    phenosv_escaped="$(escape_groovy_single_quote "$phenosv_path")"

    tmp_file="$(make_tmpfile)"
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
        echo "Error: failed to update annovar_host_path/phenosv_host_path in nextflow.config" >&2
        exit 1
    fi
    mv "$tmp_file" "$config_file"

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
        echo "Error: failed to update profiles.standard in nextflow.config" >&2
        exit 1
    fi
    mv "$tmp_file" "$config_file"
}

for arg in "$@"; do
    case "$arg" in
        light)
            version="light"
            ;;
        --non-interactive)
            non_interactive="yes"
            ;;
        --profile=*)
            selected_profile="${arg#*=}"
            ;;
        --annovar-dir=*)
            annovar_dir="${arg#*=}"
            ;;
        --phenosv-dir=*)
            phenosv_dir="${arg#*=}"
            ;;
        --annovar-bind=*)
            annovar_bind_path="${arg#*=}"
            ;;
        --phenosv-bind=*)
            phenosv_bind_path="${arg#*=}"
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "Unknown option: $arg" >&2
            usage
            exit 1
            ;;
    esac
done

current_folder="$(pwd)"

if [[ -z "$annovar_dir" ]]; then
    annovar_dir="${current_folder}/annovar"
fi
if [[ -z "$phenosv_dir" ]]; then
    phenosv_dir="${current_folder}/PhenoSV_model"
fi

if [[ "$non_interactive" == "no" ]]; then
    read -r -p "ANNOVAR directory [${annovar_dir}]: " user_annovar_dir
    read -r -p "PhenoSV download directory [${phenosv_dir}]: " user_phenosv_dir
    annovar_dir="${user_annovar_dir:-$annovar_dir}"
    phenosv_dir="${user_phenosv_dir:-$phenosv_dir}"
fi

if [[ ! -d "$annovar_dir" ]]; then
    echo "Error: ANNOVAR directory not found: $annovar_dir" >&2
    exit 1
fi
if [[ ! -f "${annovar_dir}/annotate_variation.pl" ]]; then
    echo "Error: annotate_variation.pl not found in: $annovar_dir" >&2
    exit 1
fi

mkdir -p "$phenosv_dir"

echo "Preparing ANNOVAR databases in: $annovar_dir"
pushd "$annovar_dir" >/dev/null
perl annotate_variation.pl -buildver hg38 -downdb cytoBand humandb/
perl annotate_variation.pl -buildver hg38 -downdb -webfrom annovar avsnp147 humandb/
perl annotate_variation.pl -buildver hg38 -downdb -webfrom annovar dbnsfp47a humandb/
perl annotate_variation.pl -buildver hg38 -downdb -webfrom annovar clinvar_20240917 humandb/
perl annotate_variation.pl -buildver hg38 -downdb -webfrom annovar refGene humandb/
perl annotate_variation.pl -buildver hg38 -downdb -webfrom annovar exac03 humandb/
perl annotate_variation.pl -buildver hg38 -downdb -webfrom annovar gnomad41_exome humandb/
perl annotate_variation.pl -buildver hg38 -downdb -webfrom annovar gnomad41_genome humandb/
perl annotate_variation.pl -buildver hg38 -downdb -webfrom annovar GTEx_v8_eQTL humandb/
perl annotate_variation.pl -buildver hg38 -downdb -webfrom annovar GTEx_v8_sQTL humandb/
popd >/dev/null

echo "Preparing PhenoSV resources in: $phenosv_dir"
pushd "$phenosv_dir" >/dev/null
if [[ "$version" != "light" ]]; then
    if [[ ! -f "PhenosvFile.tar" ]]; then
        echo "downloading PhenoSV files........"
        wget https://www.openbioinformatics.org/PhenoSV/PhenosvFile.tar
    fi
    echo "unzipping PhenosvFile.tar........"
    tar -xvf "PhenosvFile.tar"
    rm -f "PhenosvFile.tar"
else
    if [[ ! -f "PhenosvlightFile.tar" ]]; then
        echo "downloading PhenoSV-light files........"
        wget https://www.openbioinformatics.org/PhenoSV/PhenosvlightFile.tar
    fi
    echo "unzipping PhenosvlightFile.tar........"
    tar -xvf "PhenosvlightFile.tar"
    rm -f "PhenosvlightFile.tar"
fi

wget -O H2GKBs.zip https://github.com/WGLab/Phen2Gene/releases/download/1.1.0/H2GKBs.zip
unzip -q -o H2GKBs.zip
rm -f H2GKBs.zip

if [[ -d data ]]; then
    pushd data >/dev/null
    if [[ -f featuremaster_scu1026.csv ]]; then
        sed 's|/home/xu3/PhenoSV/data|/PhenoSV/train_data/data|g' featuremaster_scu1026.csv > features_set.csv
    fi
    if [[ -f features_set_light.csv ]]; then
        sed -i 's|/home/xu3/PhenoSV/data|/PhenoSV/train_data/data|g' features_set_light.csv
    fi
    popd >/dev/null
fi
popd >/dev/null

if [[ -z "$annovar_bind_path" ]]; then
    annovar_bind_path="$annovar_dir"
fi
if [[ -z "$phenosv_bind_path" ]]; then
    phenosv_bind_path="$phenosv_dir"
fi

if [[ "$non_interactive" == "no" && -z "$selected_profile" ]]; then
    cat <<'PROMPT'

Select your default runtime profile:
  1) standard (SLURM + Singularity)
  2) slurm_singularity (SLURM + Singularity)
  3) local_singularity (local + Singularity)
  4) local_docker (local + Docker)
PROMPT
    read -r -p "Enter choice [1-4] (default: 1): " profile_choice
    case "${profile_choice:-1}" in
        1) selected_profile="standard" ;;
        2) selected_profile="slurm_singularity" ;;
        3) selected_profile="local_singularity" ;;
        4) selected_profile="local_docker" ;;
        *) echo "Invalid profile choice: ${profile_choice}" >&2; exit 1 ;;
    esac
fi

if [[ -z "$selected_profile" ]]; then
    selected_profile="standard"
fi

if [[ "$non_interactive" == "no" ]]; then
    read -r -p "ANNOVAR bind host path [${annovar_bind_path}]: " user_annovar_bind
    read -r -p "PhenoSV bind host path [${phenosv_bind_path}]: " user_phenosv_bind
    annovar_bind_path="${user_annovar_bind:-$annovar_bind_path}"
    phenosv_bind_path="${user_phenosv_bind:-$phenosv_bind_path}"
fi

case "$selected_profile" in
    standard|slurm_singularity|local_singularity|local_docker) ;;
    *)
        echo "Invalid profile: $selected_profile" >&2
        exit 1
        ;;
esac

update_nextflow_config "$selected_profile" "$annovar_bind_path" "$phenosv_bind_path"

echo
echo "Updated config: ${current_folder}/nextflow.config"
echo "  profiles.standard backend = ${selected_profile}"
echo "  annovar_host_path= ${annovar_bind_path}"
echo "  phenosv_host_path= ${phenosv_bind_path}"
echo
echo "PipeVar will now use these defaults directly from nextflow.config (via profiles.standard)."
