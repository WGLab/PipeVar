#!/bin/bash
set -euo pipefail

version=""
non_interactive="no"
selected_profile=""
annovar_dir=""
phenosv_dir=""
annovar_bind_path=""
phenosv_bind_path=""
phenotype_image_dir=""
pull_phenotype_images="yes"

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
  --phenotype-image-dir=<path>
                            SIF destination for Singularity profiles
  --skip-phenotype-images   Do not pull PhenoGPT2 or PhenoTagger images
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
        --phenotype-image-dir=*)
            phenotype_image_dir="${arg#*=}"
            ;;
        --skip-phenotype-images)
            pull_phenotype_images="no"
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

current_folder="$(pwd -P)"
helper_directory="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"

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
# Store absolute paths so container mounts do not depend on the task directory.
annovar_dir="$(cd "$annovar_dir" && pwd -P)"
phenosv_dir="$(cd "$phenosv_dir" && pwd -P)"

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

if [[ "$selected_profile" != "local_docker" && -z "$phenotype_image_dir" ]]; then
    phenotype_image_dir="${current_folder}/containers"
fi

if [[ "$non_interactive" == "no" ]]; then
    read -r -p "ANNOVAR bind host path [${annovar_bind_path}]: " user_annovar_bind
    read -r -p "PhenoSV bind host path [${phenosv_bind_path}]: " user_phenosv_bind
    annovar_bind_path="${user_annovar_bind:-$annovar_bind_path}"
    phenosv_bind_path="${user_phenosv_bind:-$phenosv_bind_path}"
    if [[ "$selected_profile" != "local_docker" && "$pull_phenotype_images" == "yes" ]]; then
        read -r -p "Phenotype SIF directory [${phenotype_image_dir}]: " user_phenotype_image_dir
        phenotype_image_dir="${user_phenotype_image_dir:-$phenotype_image_dir}"
    fi
fi

case "$selected_profile" in
    standard|slurm_singularity|local_singularity|local_docker) ;;
    *)
        echo "Invalid profile: $selected_profile" >&2
        exit 1
        ;;
esac

# The standalone helpers own configuration updates for setup and later changes.
bash "${helper_directory}/update_bind_paths.sh" "$annovar_bind_path" "$phenosv_bind_path"
bash "${helper_directory}/update_profile.sh" "$selected_profile"

if [[ "$pull_phenotype_images" == "yes" ]]; then
    phenotype_image_args=(
        "--profile=${selected_profile}"
        "--config=${current_folder}/nextflow.config"
    )
    if [[ "$selected_profile" != "local_docker" ]]; then
        phenotype_image_args+=("--image-dir=${phenotype_image_dir}")
    fi
    bash "${helper_directory}/scripts/prepare_phenotype_images.sh" "${phenotype_image_args[@]}"
else
    echo "Skipping PhenoGPT2 and PhenoTagger image pulls."
fi

echo
echo "Updated config: ${current_folder}/nextflow.config"
echo "  profiles.standard backend = ${selected_profile}"
if [[ "$pull_phenotype_images" == "yes" && "$selected_profile" != "local_docker" ]]; then
    echo "  phenotype SIF directory   = ${phenotype_image_dir}"
fi
echo
echo "PipeVar will now use these defaults directly from nextflow.config (via profiles.standard)."
