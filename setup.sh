#!/bin/bash
set -euo pipefail

version=""
non_interactive="no"
selected_profile=""
annovar_dir=""
phenosv_dir=""
annovar_bind_path=""
phenosv_bind_path=""

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

cat > .pipevar.user.config <<EOF
/*
 * Auto-generated by setup.sh
 * You can edit this file later to change default runtime behavior.
 */
manifest {
    defaultProfile = '${selected_profile}'
}
params {
    annovar_host_path = '${annovar_bind_path}'
    phenosv_host_path = '${phenosv_bind_path}'
}
EOF

echo
echo "Wrote local user config: ${current_folder}/.pipevar.user.config"
echo "  defaultProfile   = ${selected_profile}"
echo "  annovar_host_path= ${annovar_bind_path}"
echo "  phenosv_host_path= ${phenosv_bind_path}"
echo
echo "You can now run PipeVar without repeatedly passing bind paths or -profile."
