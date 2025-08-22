# PipeVar

PipeVar is a pathogenic variant prioritization workflow for undiagnosed, rare diseases. It utilizes various tools developed from WGLab and other softwares to call structural variants, singule-nucleotide variants, indels and repeat expansions, and prioritize potential pathogenic variants and existing pathogenic variants as well.

PipeVar is implenmetned in Nextflow, and can be ran using Docker or Singularity. We are in the development of utilizing Conda for running, but for the best consistency, use either Docker or Singularity for running PipeVar.


# Usage

nextflow run main.nf --bam <FILE>
nextflow run main.nf --vcf <FILE> 


REQUIRED PARAMETERS:
  --bam <FILE> Path to input BAM file. Cannot be used with VCF option. Must be full path.
  --vcf <FILE> Path to input VCF file. Cannot be used with BAM option. Must be full path.
  --ref_fa <FILE> Reference genome in FASTA format. Must be full path.
  --out_prefix <STRING> Prefix for output files
  --note <FILE> Clinical note text file, in a format of VCF. used for HPO term extraction. Only neded if HPO terms are not available.
  --hpo <FILE> HPO ID file; note file can be used instead.


# Update that needs to be done.
