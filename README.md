# PipeVar

PipeVar is a pathogenic variant prioritization workflow for undiagnosed, rare diseases. It utilizes various tools developed from WGLab and other softwares to call structural variants, singule-nucleotide variants, indels and repeat expansions, and prioritize potential pathogenic variants and existing pathogenic variants as well.

PipeVar is implenmetned in Nextflow, and can be ran using Docker or Singularity. We are in the development of utilizing Conda for running, but for the best consistency, use either Docker or Singularity for running PipeVar. We currently only have support for Slurm, but are working on other cluster system as well.


graph TD
    subgraph Inputs
        A[Input Files: BAM/VCF, HPO, Ref Genome]
    end

    A --> B{Input Data Type?};
    B -- VCF --> C[VCF-Based Analysis];
    B -- BAM --> D[BAM-Based Analysis];

    subgraph VCF-Based Analysis
        C --> C1{Mode?};
        C1 -- --mode sv --> C2[ANNOVAR_SV];
        C2 --> C3[SURVIVOR];
        C3 --> C4[PhenoSV];
        C4 --> C_OUT[Prioritized SVs];

        C1 -- --mode snp --> C5[ANNOVAR];
        C_HPO[HPO File] --> C6[Phen2gene];
        C5 --> C7[RankVar];
        C6 --> C7;
        C5 --> C8[Rankscore_analysis]
        C6 --> C8
        C7 --> C_OUT2[Prioritized SNVs];
        C8 --> C_OUT2
    end

    subgraph BAM-Based Analysis
        D --> D1{Read Type?};
        D1 -- Long-Read (Default) --> LR[Long-Read Pipelines];
        D1 -- --type short --> SR[Short-Read Pipelines];

        subgraph Long-Read Pipelines
            LR --> LR1{Analysis Mode?};

            LR1 -- Full Analysis (Default) --> LR_FULL;
            subgraph Full Long-Read Analysis
                LR_FULL -- SNV Analysis --> LR_SNV_CALLER{--light yes?};
                LR_SNV_CALLER -- No (Default) --> CLAIR3[clair3];
                LR_SNV_CALLER -- Yes --> NANOCALLER[nanocaller];
                CLAIR3 --> ANNOVAR1[ANNOVAR];
                NANOCALLER --> ANNOVAR1;
                HPO1[HPO File] --> PHEN2GENE1[Phen2gene];
                ANNOVAR1 --> RANKVAR1[RankVar];
                PHEN2GENE1 --> RANKVAR1;
                ANNOVAR1 --> RANKSCORE1[Rankscore_analysis]
                PHEN2GENE1 --> RANKSCORE1
                RANKVAR1 --> LROUT1[Prioritized SNVs];
                RANKSCORE1 --> LROUT1
                
                LR_FULL -- SV Analysis --> CUTESV[cuteSV];
                LR_FULL -- SV Analysis --> SNIFFLES[sniffles];
                CUTESV --> TRUVARI[truvari];
                SNIFFLES --> TRUVARI;
                TRUVARI --> SURVIVOR1[SURVIVOR];
                SURVIVOR1 --> PHENOSV1[PhenoSV];
                PHENOSV1 --> LROUT2[Prioritized SVs];

                LR_FULL -- STR Analysis --> NANOREPEAT[NanoRepeat];
                NANOREPEAT --> LROUT3[STR Calls];
            end

            LR1 -- --mode sv --> LR_SV_ONLY[SV Analysis Only];
            subgraph SV Only
                LR_SV_ONLY --> CUTESV2[cuteSV];
                LR_SV_ONLY --> SNIFFLES2[sniffles];
                CUTESV2 --> TRUVARI2[truvari];
                SNIFFLES2 --> TRUVARI2;
                TRUVARI2 --> SURVIVOR2[SURVIVOR];
                SURVIVOR2 --> PHENOSV2[PhenoSV];
                PHENOSV2 --> LROUT4[Prioritized SVs];
            end

            LR1 -- --mode snp --> LR_SNP_ONLY[SNP Analysis Only];
            subgraph SNP Only
                LR_SNP_ONLY --> LR_SNV_CALLER2{--light yes?};
                LR_SNV_CALLER2 -- No (Default) --> CLAIR3_2[clair3];
                LR_SNV_CALLER2 -- Yes --> NANOCALLER_2[nanocaller];
                CLAIR3_2 --> ANNOVAR2[ANNOVAR];
                NANOCALLER_2 --> ANNOVAR2;
                HPO2[HPO File] --> PHEN2GENE2[Phen2gene];
                ANNOVAR2 --> RANKVAR2[RankVar];
                PHEN2GENE2 --> RANKVAR2;
                ANNOVAR2 --> RANKSCORE2[Rankscore_analysis]
                PHEN2GENE2 --> RANKSCORE2
                RANKVAR2 --> LROUT5[Prioritized SNVs];
                RANKSCORE2 --> LROUT5
            end
        end

        subgraph Short-Read Pipelines
            SR --> SR1{Lightweight Mode?};
            SR1 -- No (Default) --> SR_DEFAULT[Full Short-Read Analysis];
            subgraph Full Short-Read Analysis
                SR_DEFAULT -- SNV Analysis --> DEEPVARIANT[deepvariant];
                DEEPVARIANT --> ANNOVAR3[ANNOVAR];
                HPO3[HPO File] --> PHEN2GENE3[Phen2gene];
                ANNOVAR3 --> RANKVAR3[RankVar];
                PHEN2GENE3 --> RANKVAR3;
                RANKVAR3 --> SROUT1[Prioritized SNVs];
                
                SR_DEFAULT -- SV Analysis --> MANTA[Manta];
                MANTA --> SURVIVOR3[SURVIVOR];
                SURVIVOR3 --> PHENOSV3[PhenoSV];
                PHENOSV3 --> SROUT2[Prioritized SVs];

                SR_DEFAULT -- STR Analysis --> EXPHUNTER[ExpansionHunter];
                EXPHUNTER --> SROUT3[STR Calls];
            end

            SR1 -- --light yes --> SR_LIGHT[Light Short-Read Analysis];
             subgraph Light Short-Read Analysis
                SR_LIGHT -- SNV Analysis --> HAPLOTYPECALLER[haplotypecaller];
                HAPLOTYPECALLER --> ANNOVAR4[ANNOVAR];
                HPO4[HPO File] --> PHEN2GENE4[Phen2gene];
                ANNOVAR4 --> RANKVAR4[RankVar];
                PHEN2GENE4 --> RANKVAR4;
                RANKVAR4 --> SROUT4[Prioritized SNVs];

                SR_LIGHT -- SV Analysis --> MANTA2[Manta];
                MANTA2 --> SURVIVOR4[SURVIVOR];
                SURVIVOR4 --> PHENOSV4[PhenoSV];
                PHENOSV4 --> SROUT5[Prioritized SVs];

                SR_LIGHT -- STR Analysis --> EXPHUNTER2[ExpansionHunter];
                EXPHUNTER2 --> SROUT6[STR Calls];
            end
        end
    end
# Requirements

PipeVar requires either Docker or Singularity to run. If your system do not have Singularity installed as a module, you can try to install Singularity using conda with


  conda create -n singularity singularity

  or install singualrity in your conda environment by

  conda install conda-forge::singularity


# Set up

To use PipeVar, there is a set up stage required for two softwares. The first one is ANNOVAR.

To download ANNOVAR,


You will need to download PhenoSV model files. You can download the light version, or download the full version. If you download the light version, the pipeline can be only ran with --light yes option.

After downloading PhenoSV model file, run this setup.sh file in your folder. It will modify nextflow.config to contain your required phenosv and download necessary datasets to run ANNOVAR.


# Usage

nextflow run main.nf --bam <FILE>
nextflow run main.nf --vcf <FILE> 


  REQUIRED PARAMETERS:
  
  --bam <FILE> Path to input BAM file. Cannot be used with VCF option. Must be full path. Requires .bai index file.
  
  --vcf <FILE> Path to input VCF file. Cannot be used with BAM option. Must be full path.
  
  --ref_fa <FILE> Reference genome in FASTA format. Must be full path.
  
  --out_prefix <STRING> Prefix for output files
  
  --note <FILE> Clinical note text file, in a format of VCF. used for HPO term extraction. Only neded if HPO terms are not available.
  
  --hpo <FILE> HPO ID file; note file can be used instead.


  # Example hpo.txt
    HP:0001250
    HP:0000750
    HP:0001257

OPTIONAL PARAMETERS:
    --input_directory <DIR>   Directory containing input files
    --output_directory <DIR>  Path to output directory (default: current dir)
    --mode <sv|snp>           Variant type to analyze (required with --vcf or --bam)
    --type <short|long>       Input data type: short or long reads (required with --bam)
    --light <yes|no>          Use lightweight PhenoSV model and NanoCaller (faster, lower memory)
    --gq <INT>                Minimum genotype quality [default: 20]
    --ad <INT>                Minimum allelic depth [default: 15]
    --gnomad <FLOAT>          Max gnomAD allele frequency [default: 0.0001]
    --help                    Print this help message and exit

EXAMPLES:

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

NOTES:
    - At least one of `--hpo` or `--note` must be provided.
    - If `--note` is used, `--hpo` is auto-generated via phenotagger.
    - `--mode` must be specified for VCF input, and helps direct SNV vs SV flow.
    - `--type` is required for BAM input to specify sequencing technology.
    - All file paths must be absolute or relative to `--input_directory`.
    - `--light yes` uses faster, resource-friendly software such as haplotypecaller, NanoCaller and PhenoSV-light.

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

Add LongPhase process, and ACMG Guideline, and PhenoGPT2 for note direction.



