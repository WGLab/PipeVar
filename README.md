# PipeVar

PipeVar is a pathogenic variant prioritization workflow for undiagnosed, rare diseases. It utilizes various tools developed from WGLab and other softwares to call structural variants, singule-nucleotide variants, indels and repeat expansions, and prioritize potential pathogenic variants and existing pathogenic variants as well.

PipeVar is implenmetned in Nextflow, and can be ran using Docker or Singularity. We are in the development of utilizing Conda for running, but for the best consistency, use either Docker or Singularity for running PipeVar. We currently only have support for Slurm, but are working on other cluster system as well.

<img width="3840" height="2182" alt="Untitled diagram _ Mermaid Chart-2025-08-29-153706" src="https://github.com/user-attachments/assets/d0c1b7d2-0dc8-49ef-8da5-b86b3b51aa01" />


# Requirements

PipeVar requires either Docker or Singularity to run. If your system do not have Singularity installed as a module, you can try to install Singularity using conda with

```
conda create -n singularity singularity
```
  or install singualrity in your conda environment by
```
conda install conda-forge::singularity
```

Currently, PipeVar is only tested to run in SLURM environment. It can be ran in other environment by changing executor paramemter in nextflow.config, but it has not been tested yet. We plan to test PipeVar in different environments.

# Set up

First, download the git repo using:

```
#Download PipeVar
git clone https://github.com/WGLab/PipeVar.git
```

To use PipeVar, there is a set up stage required for two software, ANNOVAR and PhenoSV.

To download ANNOVAR, go to this link https://www.openbioinformatics.org/annovar/annovar_download_form.php , and follow the instruction in the format. Make sure to have ANNOVAR in the PipeVar folder for setup script to run properly.

Once ANNOVAR is downloaded, run setup.sh in the directory you cloned PipeVar to download the necessary files for ANNOVAR and PhenoSV using
```
#Install full PhenoSV.
./setup.sh 
```
or
```
#Install lighter version of PhenoSV.
./setup.sh light
```
If you want to download a lighter version of PhenoSV. PhenoSV light is about ~50GB in size, while full PhenoSV is about ~150GB in size. The setup script will also modify the nextflow.config to annovar/PhenoSV directory location to necessary location to run the pipeline.

# Usage
```
#BAM version.
nextflow run main.nf --bam <FILE> --ref_fa <FILE> --out_prefix <FOLDER> --note <FILE> or hpo <FILE>

#VCF version. Requires either SV or SNV option.
nextflow run main.nf --vcf <FILE> --mode sv or snv --ref_fa <FILE> --out_prefix <FOLDER> --note <FILE> or hpo <FILE> 
```

```
REQUIRED PARAMETERS:
  
  --bam <FILE> Path to input BAM file. Cannot be used with VCF option. Must be full path. Requires .bai index file.
  
  --vcf <FILE> Path to input VCF file. Cannot be used with BAM option. Must be full path.
  
  --ref_fa <FILE> Reference genome in FASTA format. Must be full path.
  
  --out_prefix <STRING> Prefix for output files
  
  --note <FILE> Clinical note text file, in a format of VCF. used for HPO term extraction. Only neded if HPO terms are not available.
  
  --hpo <FILE> HPO ID file; note file can be used instead.

  --mode <sv|snv>. Option run either SV mode or SNV mode. Required for VCF mode. Optional for BAM.

OPTIONAL PARAMETERS:
    --output_directory <DIR>  Path to output directory (default: current directory)
    --type <ont|pacbio|short>       Input data type: short/long reads(either Pac-Bio or ONT) (default is ONT).
    --light <yes|no>          Use lightweight PhenoSV model, NanoCaller (faster, lower memory, but with lower accuracy)
    --gq <INT>                Minimum genotype quality [default: 20] used for filtering for RankVar and RankScore analysis.
    --ad <INT>                Minimum allelic depth [default: 15] used for filtering for RankVar and RankScore analysis.
    --gnomad <FLOAT>          Max gnomAD allele frequency [default: 0.0001] used for filtering for RankVar and RankScore analysis.
    --help                    Print this help message and exit
```

EXAMPLES:

```
    1. Long-read full pipeline (SV + SNP + STR):
        nextflow run main.nf \
          --bam /data/sample.bam \
          --ref_fa /refs/hg38.fa \
          --out_prefix patient1 \
          --hpo /data/hpo.txt \
          --type long

    2. Short-read full pipeline:
        nextflow run main.nf \
          --bam /data/sample.bam \
          --ref_fa /refs/hg38.fa \
          --out_prefix patient1 \
          --hpo /data/hpo.txt \
          --type short

    3. Short-read with lightweight model:
        nextflow run main.nf \
          --bam /data/sample.bam \
          --ref_fa /refs/hg38.fa \
          --out_prefix patient1 \
          --hpo /data/hpo.txt \
          --type short \
          --light yes

    4. Variant re-annotation using VCF (SV mode):
        nextflow run main.nf \
          --vcf /data/sample.vcf \
          --ref_fa /refs/hg38.fa \
          --out_prefix patient_sv \
          --hpo /data/hpo.txt \
          --mode sv

    5. Auto-extract HPO from clinical notes:
        nextflow run main.nf \
          --bam /data/sample.bam \
          --ref_fa /refs/hg38.fa \
          --out_prefix patient1 \
          --note /data/note.txt \
          --type ont
```

Ideally, the job would be the best ran using the job submission since variant calling process can take long depending on your job. Following is aa simple example for SLURM job submission for default setting.

ADD EXAMPLE HERE!

NOTES:

    - At least one of `--hpo` or `--note` must be provided.
    
    - If `--note` is used, `--hpo` is auto-generated via phenotagger.
    
    - `--mode` must be specified for VCF input, and helps direct SNV vs SV flow.
    
    - `--type` is required for BAM input to specify sequencing technology.
    
    - All file paths must be absolute or relative to `--input_directory`.
    
    - `--light yes` uses faster, resource-friendly software such as haplotypecaller, NanoCaller and PhenoSV-light.

# Parameter details


***Reference genome***
```
--ref_fa
```

Used to specifcy the reference genome FASTA file. It must be indexed.

FASTA file can be indexed using samtools with following command if needed:

```
samtools faidx file.fa
``` 


***BAM file***
```
--bam
```

BAM file needed to run the workflow. Must be indexed and sorted. SAM/CRAM are not accepted for now, but will be in future. If bam files is not indexed or sorted, use following example command

```
#Sorting bam
samtools sort -o your.sorted.bam your.bam

#Index bam
samtools index -b your.bam (or your.sorted.bam if you sorted).
```

***VCF file***
```
--vcf
```
VCF file can be used if you want to run PipeVar on SV or SNV. Does not require indexing or sorting.

***Out prefix***

```
--out_prefix
```
String for the prefix for the output and processing. Makes sure to not include any special characters or space in between.

***Medical Note***
```
--note
```
Medical notes that can be in csv, tsv or txt format. Make sure to input one medical not for one patient.

***HPO file***
```
--hpo

#Example
HP:0031647
HP:0031647

```
Text file with list of HPO terms. Make sure each HPO term is in new line, and only contain HPO terms. 

***Mode selection***

```
--mode
```

Option for either SV or SNV. Required for VCF option, but optional for BAM option. Only 'sv' or 'snv' is accepted.

***Sequencing type***

```
--type
```

Option for sequencing type. The default is 'ONT' for now, and has 'short' for NGS and 'pacbio' for Pac-Bio option. For long-read sequencing, they are used in repeat and SV process.

***Output direcrtory***
```
--output_directory
```
Option for output directory location. Must provide full path. The default is current directory where the script is submitted.


***Gnomad frequency***
```
---gnomad
```
Option for gnomad allele frequency filtering for SNV Rankscore and SV score. The default is 0.0001, but may be lowered for autosomal recessive variant priortization. SV filtering in progress.


# Softwares used

PIPELINE MODULES:

    SNV CALLING
      - clair3         : Deep learning SNP caller (long-read)
      - nanocaller     : Lightweight long-read SNP caller
      - haplotypecaller: Short-read SNP caller (GATK)
      - deepvariant    : Deep learning short-read SNP caller

    SV CALLING
      - cuteSV         : Long-read SV caller
      - sniffles       : Long-read SV caller
      - Manta          : Short-read SV caller
      - truvari        : SV comparison/benchmarking
      - SURVIVOR       : SV merging

    STR DETECTION
      - NanoRepeat     : Long-read STR caller
      - ExpansionHunter: Short-read STR detection

    PHENOTYPING
      - Phen2gene      : HPO-to-gene mapping
      - phenotagger    : NLP-based clinical note to HPO term conversion
      - PhenoGpt2      : To be implemented

    VARIANT RANKING
      - ANNOVAR        : SNV/SV annotation
      - RankVar        : Final SNV ranking
      - Rankscore_analysis: Additional ranking analysis


# Output

All the output will be stored in output folder, or the launch folder. The list of outputs are as followed:




# Update that needs to be done.

Add LongPhase process, and ACMG Guideline, and PhenoGPT2 for note direction. Clean up the output results as well.

We are also working on version so reference genome may not be previously downloaded, but can be ran with hg38 or grch38 option.

SAM/CRAM support will be included in future as well.
