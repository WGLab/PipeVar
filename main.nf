


def helpMessage = """
──────────────────────────────────────────────────────────────
🧬 PipeVar: Variant Prioritization Pipeline for Rare Diseases
──────────────────────────────────────────────────────────────

DESCRIPTION:
    PipeVar is a modular Nextflow pipeline designed to prioritize 
    disease-causal variants by integrating SNVs, SVs, STRs, and 
    phenotype (HPO) information from long-read or short-read 
    whole-genome sequencing data.

USAGE:
    nextflow run pipevar.nf [options]

REQUIRED PARAMETERS:
    --bam <FILE>              Path to input BAM file (for SNV/STR calling)
    --vcf <FILE>              Path to input VCF file (optional if BAM provided)
    --out_prefix <STRING>     Prefix for all output files and directories
    --ref_fa <FILE>           Path to reference genome FASTA file
    --hpo <FILE>              Path to HPO term list file (txt, comma- or line-separated)
    --input_directory <DIR>   Directory containing input files and context

OPTIONAL PARAMETERS:
    --gq <INT>                Minimum genotype quality (GQ) filter [default: 20]
    --ad <INT>                Minimum allelic depth (AD) filter [default: 15]
    --gnomad <FLOAT>          Maximum gnomAD allele frequency [default: 0.0001]
    --light <yes|no>          Use lightweight PhenoSV model [default: null]
    --mode <STRING>           Mode flag (reserved for future use)
    --type <STRING>           Type flag (reserved for future use)
    --help                    Print this help message and exit

EXAMPLES:
    Run full pipeline with a BAM file:
        nextflow run pipevar.nf \
          --bam /path/to/sample.bam \
          --out_prefix patient1_results \
          --ref_fa /refs/hg38.fa \
          --hpo /path/to/hpo_terms.txt \
          --input_directory /project/input

    Run pipeline with light model and custom filters:
        nextflow run pipevar.nf \
          --bam /data/sample.bam \
          --out_prefix sample_light \
          --ref_fa /refs/hg38.fa \
          --hpo /data/hpo.txt \
          --input_directory /data \
          --light yes \
          --gq 25 --ad 20 --gnomad 0.001

PIPELINE MODULES INCLUDED:
    - cuteSV         : Structural variant detection
    - sniffles       : SV calling with support for long reads
    - nanocaller     : Variant calling using long reads
    - clair3         : Deep learning-based variant calling
    - SURVIVOR       : SV merging and comparison
    - PhenoSV        : Phenotype-prioritized SV ranking
    - ANNOVAR        : Variant annotation
    - Phen2gene      : Gene prioritization from HPO terms
    - RankVar        : Final variant ranking
    - NanoRepeat     : STR detection from ONT reads

RESOURCE USAGE (SLURM):
    Each process is configured with SLURM:
      - nanocaller   : 4 CPUs, 32 GB
      - PhenoSV      : 2 CPUs, 16 GB, 48h
      - clair3       : 4 CPUs, 32 GB, 24h
      - NanoRepeat   : 4 CPUs, 32 GB
      - ANNOVAR      : 2 CPUs, 16 GB
      - Phen2gene    : 2 CPUs, 8 GB

NOTES:
    - All paths should be absolute.
    - Input files must be placed within or referenced relative to --input_directory.
    - The `--light yes` option activates a faster, less resource-intensive PhenoSV model.

──────────────────────────────────────────────────────────────


"""
if (params.help) {
	println(helpMessage)
	exit(0)
}


	


include { cuteSV } from './modules/cutesv/'
include { sniffles } from './modules/sniffles/'
include { nanocaller } from './modules/nanocaller/'
include { clair3 } from './modules/clair3/'
include { SURVIVOR } from './modules/survivor/'
include { PhenoSV } from './modules/phenosv/'
include { ANNOVAR } from './modules/annovar/'
include { Phen2gene } from './modules/phen2gene/'
include { RankVar } from './modules/rankvar/'
include { NanoRepeat } from './modules/nanorepeat/'
include { truvari } from './modules/truvari/'
include { Manta } from './modules/manta/'
include { haplotypecaller } from './modules/haplotypecaller/'
include { deepvariant } from './modules/deepvariant/'
include { ANNOVAR_SV } from './modules/annovar_sv/'
include { ExpansionHunter } from './modules/expansion_hunter/'
include { Rankscore_analysis } from './modules/rankscore_analysis/'
include { phenotagger } from './modules/phenotagger/'
include { longphase } from './modules/longphase/'



workflow {
	if ( params.bam ) {
		bam=Channel.fromPath(params.bam, checkIfExists: true)
		input_directory=bam.map { it.parent.toString() }
	}
	else if ( params.vcf ) { 
		vcf=Channel.fromPath(params.vcf, checkIfExists: true)
		input_directory=vcf.map { it.parent.toString() }
	}
	if ( params.ref_fa ) {
		ref_fa=Channel.fromPath(params.ref_fa, checkIfExists: true)
		ref_fa_directory=ref_fa.map {it.parent.toString() }
	}

	if ( params.note ) {
		note=Channel.fromPath(params.note, checkIfExists: true)
                note_directory=note.map {it.parent.toString() }
	}
	if ( params.hpo ) { 
		hpo=Channel.fromPath(params.hpo, checkIfExists: true)
		hpo_directory=hpo.map {it.parent.toString() }
	}
	output_directory_check=file(params.output_directory)
	if ( !output_directory_check.exists() ) {
                output_directory_check.mkdirs()
        }
	output_directory=Channel.fromPath(params.output_directory)
	output_directory_full=Channel.value(params.output_directory)		
	out_prefix=Channel.value(params.out_prefix)
	type=Channel.value(params.type)
	gnomad=Channel.value(params.gnomad)
	gq=Channel.value(params.gq)
	ad=Channel.value(params.ad)
	if ( params.note && !params.hpo ) {
		phenotagger(note,note_directory,output_directory,output_directory_full,out_prefix)
		hpo=phenotagger.out
		hpo_directory=Channel.value(params.output_directory)
	}
	if ( params.vcf && params.mode == 'sv' ) {
		ANNOVAR_SV(vcf,out_prefix,input_directory,output_directory,output_directory_full)
		SURVIVOR(ANNOVAR_SV.out,out_prefix,output_directory,output_directory_full)
		PhenoSV(SURVIVOR.out,out_prefix,hpo,hpo_directory,output_directory,output_directory_full)
	}
	else if ( params.vcf && params.mode == 'snp' ) {
		ANNOVAR(vcf,input_directory,out_prefix,output_directory,output_directory_full)
                Phen2gene(hpo,hpo_directory,out_prefix,output_directory,output_directory_full)
                RankVar(ANNOVAR.out,Phen2gene.out,hpo,hpo_directory,out_prefix,gnomad,gq,ad,output_directory,output_directory_full)
		Rankscore_analysis(ANNOVAR.Out,Phen2gene.out,out_prefix,output_directory,output_directory_full)
	}
	else if ( params.bam && params.mode == 'sv' ) {
		cuteSV(bam,input_directory,out_prefix,ref_fa,ref_fa_directory,output_directory,output_directory_full)
                sniffles(bam,input_directory,out_prefix,ref_fa,ref_fa_directory,output_directory,output_directory_full)
                truvari(cuteSV.out,sniffles.out,ref_fa,ref_fa_directory,out_prefix,output_directory,output_directory_full)
                SURVIVOR(truvari.out,out_prefix,output_directory,output_directory_full)
                PhenoSV(SURVIVOR.out,out_prefix,hpo,hpo_directory,output_directory,output_directory_full)
	}
        else if ( params.bam && params.mode == 'snp' && params.light == 'yes' ) {
                nanocaller(bam,input_directory,out_prefix,ref_fa,ref_fa_directory,output_directory,output_directory_full)
                ANNOVAR(nanocaller.out,output_directory_full,out_prefix,output_directory,output_directory_full)
                Phen2gene(hpo,hpo_directory,out_prefix,output_directory,output_directory_full)
                RankVar(ANNOVAR.out,Phen2gene.out,hpo,hpo_directory,out_prefix,gnomad,gq,ad,output_directory,output_directory_full)
                Rankscore_analysis(ANNOVAR.Out,Phen2gene.out,out_prefix,output_directory,output_directory_full)

        }
	else if ( params.bam && params.mode == 'snp' ) {
		clair3(bam,input_directory,out_prefix,ref_fa,ref_fa_directory,output_directory,output_directory_full)
		ANNOVAR(clair3.out,output_directory_full,out_prefix,output_directory,output_directory_full)
                Phen2gene(hpo,hpo_directory,out_prefix,output_directory,output_directory_full)
                RankVar(ANNOVAR.out,Phen2gene.out,hpo,hpo_directory,out_prefix,gnomad,gq,ad,output_directory,output_directory_full)
                Rankscore_analysis(ANNOVAR.Out,Phen2gene.out,out_prefix,output_directory,output_directory_full)	
	}
	else if (params.bam && params.type == 'short' && params.light == 'yes') {
                haplotypecaller(bam,input_directory,out_prefix,ref_fa,ref_fa_directory,output_directory,output_directory_full)
                ANNOVAR(haplotypecaller.out,output_directory_full,out_prefix,output_directory,output_directory_full)
                Phen2gene(hpo,hpo_directory,out_prefix,output_directory,output_directory_full)
                RankVar(ANNOVAR.out,Phen2gene.out,hpo,hpo_directory,out_prefix,gnomad,gq,ad,output_directory,output_directory_full)
		Manta(bam,input_directory,out_prefix,ref_fa,ref_fa_directory,output_directory,output_directory_full)
                SURVIVOR(Manta.out,out_prefix,output_directory,output_directory_full)
                PhenoSV(SURVIVOR.out,out_prefix,hpo,hpo_directory,output_directory,output_directory_full)
		ExpansionHunter(bam,input_directory,out_prefix,ref_fa,ref_fa_directory,output_directory,output_directory_full)
        }
        else if (params.bam && params.type == 'short' ) {
                deepvariant(bam,input_directory,out_prefix,ref_fa,ref_fa_directory,output_directory,output_directory_full)
		ANNOVAR(deepvariant.out,output_directory_full,out_prefix,output_directory,output_directory_full)
                Phen2gene(hpo,hpo_directory,out_prefix,output_directory,output_directory_full)
                RankVar(ANNOVAR.out,Phen2gene.out,hpo,hpo_directory,out_prefix,gnomad,gq,ad,output_directory,output_directory_full)
                Manta(bam,input_directory,out_prefix,ref_fa,ref_fa_directory,output_directory,output_directory_full)
                SURVIVOR(Manta.out,out_prefix,output_directory,output_directory_full)
                PhenoSV(SURVIVOR.out,out_prefix,hpo,hpo_directory,output_directory,output_directory_full)
                ExpansionHunter(bam,input_directory,out_prefix,ref_fa,ref_fa_directory,output_directory,output_directory_full)

        }
	else if (params.bam && params.type != 'short' && params.light == 'yes') {
		nanocaller(bam,input_directory,out_prefix,ref_fa,ref_fa_directory,output_directory,output_directory_full)
                ANNOVAR(nanocaller.out,output_directory_full,out_prefix,output_directory,output_directory_full)
                Phen2gene(hpo,hpo_directory,out_prefix,output_directory,output_directory_full)
                RankVar(ANNOVAR.out,Phen2gene.out,hpo,hpo_directory,out_prefix,gnomad,gq,ad,output_directory,output_directory_full)
                Rankscore_analysis(ANNOVAR.Out,Phen2gene.out,out_prefix,output_directory,output_directory_full)
		cuteSV(bam,input_directory,out_prefix,ref_fa,ref_fa_directory,output_directory,output_directory_full)
                sniffles(bam,input_directory,out_prefix,ref_fa,ref_fa_directory,output_directory,output_directory_full)
                truvari(cuteSV.out,sniffles.out,ref_fa,ref_fa_directory,out_prefix,output_directory,output_directory_full)
                SURVIVOR(truvari.out,out_prefix,output_directory,output_directory_full)
                PhenoSV(SURVIVOR.out,out_prefix,hpo,hpo_directory,output_directory,output_directory_full)
                NanoRepeat(bam,input_directory,out_prefix,ref_fa,ref_fa_directory,output_directory,output_directory_full)
	}
	else {
		clair3(bam,input_directory,out_prefix,ref_fa,ref_fa_directory,output_directory,output_directory_full)
                ANNOVAR(clair3.out,output_directory_full,out_prefix,output_directory,output_directory_full)
		Phen2gene(hpo,hpo_directory,out_prefix,output_directory,output_directory_full)
                RankVar(ANNOVAR.out,Phen2gene.out,hpo,hpo_directory,out_prefix,gnomad,gq,ad,output_directory,output_directory_full)
                rankscore_result=Rankscore_analysis(ANNOVAR.Out,Phen2gene.out,out_prefix,output_directory,output_directory_full)
                cuteSV(bam,input_directory,out_prefix,ref_fa,ref_fa_directory,output_directory,output_directory_full)
                sniffles(bam,input_directory,out_prefix,ref_fa,ref_fa_directory,output_directory,output_directory_full)
                truvari(cuteSV.out,sniffles.out,ref_fa,ref_fa_directory,out_prefix,output_directory,output_directory_full)
                SURVIVOR(truvari.out,out_prefix,output_directory,output_directory_full)
                PhenoSV(SURVIVOR.out,out_prefix,hpo,hpo_directory,output_directory,output_directory_full)
                NanoRepeat(bam,input_directory,out_prefix,ref_fa,ref_fa_directory,output_directory,output_directory_full)
		longphase(bam,input_directory,clair3.out,cuteSV.out,output_directory,output_directory_full,PhenoSV.out,rankscore_result.rankscore,rankscore_result.clinvar,RankVar.out,out_prefix,ref_fa,ref_fa_directory)
	}
}

workflow.onComplete {
	println "Pipeline completed at: $workflow.complete"
	println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}

workflow.onError {
	println "Error: Pipeline execution stopped with the following message: ${workflow.errorMessage}"
}



