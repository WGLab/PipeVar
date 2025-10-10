#!/bin/bash


version="$1"

#After downloading annovar, download all the files necessary, and do the add annotation options

check_annovar="annovar"

# Check if the directory does NOT exist
if ! [ -d "$check_annovar" ]; then
    echo "Error: Directory '$check_annovar' not found. Aborting." >&2
    exit 1
fi

cd annovar

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

cd ..



#Download phenosv, and change the paths as well


mkdir -p PhenoSV_model

cd PhenoSV_model

if [[ "$version" != "light" ]]
then
  if ! test -f "PhenosvFile.tar"; then
  echo "downloading PhenoSV files........"
  wget https://www.openbioinformatics.org/PhenoSV/PhenosvFile.tar
  fi
  echo "unzipping PhenosvFile.tar........"
  tar -xvf "PhenosvFile.tar"
  rm "PhenosvFile.tar"
else
  if ! test -f "PhenosvlightFile.tar"; then
  echo "downloading PhenoSV-light files........"
  wget https://www.openbioinformatics.org/PhenoSV/PhenosvlightFile.tar
  fi
  echo "unzipping PhenosvlightFile.tar........"
  tar -xvf "PhenosvlightFile.tar"
  rm "PhenosvlightFile.tar"
fi

wget https://github.com/WGLab/Phen2Gene/releases/download/1.1.0/H2GKBs.zip

unzip -q H2GKBs.zip
rm H2GKBs.zip

cd data

sed 's|/home/xu3/PhenoSV/data|/PhenoSV/train_data/data|g' featuremaster_scu1026.csv > features_set.csv

sed -i 's|/home/xu3/PhenoSV/data|/PhenoSV/train_data/data|g' features_set_light.csv

cd ../..

current_folder=$(pwd)

# replace placeholder for PhenoSV_model
sed -i "s|PHENOSV_PLACEHOLDER|${current_folder}/PhenoSV_model:/PhenoSV/train_data|" nextflow.config

# replace placeholder for annovar
sed -i "s|ANNOVAR_PLACEHOLDER|${current_folder}/annovar:/annovar|" nextflow.config

