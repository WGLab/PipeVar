#!/bin/bash


version="$1"

#After downloading annovar, download all the files necessary, and do the add annotation options




#Download phenosv, and change the paths as well

#cd ..

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
