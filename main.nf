// Start here: help → input preparation → biological route. See workflow.md.
def helpMessage() {
def pipelineVersion = params.pipeline_version ?: "0.5.0"

"""
================================================================================
  PipeVar_simplified
  Version: ${pipelineVersion}
================================================================================

Rare-disease variant prioritization for short-read, long-read, and optional
mitochondrial analysis.

USAGE
  Single BAM/CRAM:
    nextflow run main.nf --bam sample.bam --ref_fa ref.fa (--note note.txt | --hpo hpo.txt) [options]

  Single VCF:
    nextflow run main.nf --vcf sample.vcf --ref_fa ref.fa --mode <snp|sv> (--note note.txt | --hpo hpo.txt) [options]

  Annotated SNV:
    nextflow run main.nf --annotated_snv yes --annovar_txt sample.hg38_multianno.txt --vcf sample.hg38_multianno.vcf (--note note.txt | --hpo hpo.txt) [options]

  Legacy CSV:
    nextflow run main.nf --input_csv samples.csv --bam true --ref_fa ref.fa [options]
    nextflow run main.nf --input_csv samples.csv --vcf true --ref_fa ref.fa --mode <snp|sv> [options]

  Unified CSV:
    nextflow run main.nf --input_csv samples.csv --ref_fa ref.fa [options]

PROFILES
  standard              SLURM + Singularity
  slurm_singularity     SLURM + Singularity
  local_singularity     Local executor + Singularity
  local_docker          Local executor + Docker

COMMON OPTIONS
  --type <ont|pacbio|short>     Sequencing type for BAM/CRAM flows
  --mode <snp|sv>               Restrict to one branch where supported
  --light <yes|no>              Use lightweight callers/models where supported
  --mito <yes|no>               Enable mitochondrial analysis for BAM/CRAM input
  --xtea <yes|no>               Enable short-read mobile-element calling
  --GPU <yes|no>                Enable shared GPU mode for DeepVariant GPU and PhenoGPT2
  --phenotype_extractor <STR>   phenotagger or phenogpt2 for clinical notes
  --phenogpt2_model_host_path   Complete versioned new_model directory (PhenoGPT2 notes)
  --phenogpt2_negation_model_host_path   Complete negation model directory
  --phenogpt2_embedding_model_host_path  Complete embedding model directory
  --phenogpt2_cache_host_path   Optional pre-created persistent cache directory
  --phenogpt2_sif <FILE>        Pre-staged PhenoGPT2 SIF (Singularity only)
  --phenotagger_sif <FILE>      Pre-staged PhenoTagger SIF (Singularity only)
  --out_prefix <STRING>         Single-sample output prefix
  --output_directory <DIR>      Publish directory
  --nanocaller_dp <NUMBER>      Missing-GQ NanoCaller depth threshold (default: 20)
  --gnomad_af_ad <FLOAT>        Dominant/XLD AF ceiling (default: 0.001)
  --gnomad_af_ar <FLOAT>        Recessive/XLR and upstream AF ceiling (default: 0.01)
  --common_sv_af_ad <FLOAT>     Dominant/XLD/unknown SV ceiling (default: 0.005)
  --common_sv_af_ar <FLOAT>     Recessive/XLR and upstream SV ceiling (default: 0.01)
  --include_clinvar_report      Include/exclude ClinVar-exclusive report entries;
                                P/LP classification is unchanged (default: yes)
  --allow_unphased_comphet      Permit unresolved AR pairs (default: no); PS is
                                optional and missing PS falls back to phased GT

NOTES
  - BAM/CRAM inputs require index files (.bai or .crai).
  - Reference FASTA index (.fai) must exist when --ref_fa is supplied.
  - Single-file mode requires one phenotype source: --note <FILE> or --hpo <FILE>.
  - PhenoGPT2 clinical-note runs require an external read-only model mount;
    HPO-only inputs do not require a PhenoGPT2 model, cache, or GPU.
  - PhenoSV receives simple DEL/DUP/INV/INS events as BED and BND/TRA
    adjacencies as BEDPE. Light mode warns that translocation accuracy is reduced.
  - Valid PS metadata is retained when present. Missing or malformed PS uses GT
    orientation; different valid phase sets remain unresolved.
  - Detailed input schemas, parameter defaults, examples, and outputs are in README.md.
================================================================================
"""
}

def validateStagedSifEngine(containerEngine, phenogpt2Sif, phenotaggerSif) {
    def supplied = [phenogpt2_sif: phenogpt2Sif, phenotagger_sif: phenotaggerSif]
        .findAll { option, value -> value != null && value.toString().trim() }
    if (!supplied) {
        return
    }
    def engine = containerEngine?.toString()?.trim()?.toLowerCase()
    if (!(engine in ['singularity', 'apptainer'])) {
        def options = supplied.keySet().collect { "--${it}" }.join(', ')
        throw new IllegalArgumentException("${options} may be used only with a Singularity profile; active container engine: ${engine ?: 'none'}.")
    }
}

include { INPUT_PREPARATION } from './subworkflows/input_preparation/main'
include { ALIGNMENT_LONG } from './subworkflows/long/main'
include { ALIGNMENT_NGS } from './subworkflows/ngs/main'
include { ALIGNMENT_LONG_MITO } from './subworkflows/mitochondrial/long'
include { ALIGNMENT_NGS_MITO } from './subworkflows/mitochondrial/ngs'
include { ALIGNMENT_VCF } from './subworkflows/vcf/main'
include { PREANNOTATED_SMALL_VARIANT } from './subworkflows/preannotated_snp/core'
include { PREANNOTATED_SV } from './subworkflows/preannotated_sv/main'
include { PREANNOTATED_ALL_NGS } from './subworkflows/preannotated_combined/main'
include { PREANNOTATED_SNV_CALLED_SV_NGS } from './subworkflows/preannotated_combined/called_sv'

workflow {
    validateStagedSifEngine(workflow.containerEngine, params.phenogpt2_sif, params.phenotagger_sif)

    if (params.help) {
        println(helpMessage())
        return
    }

    prepared = INPUT_PREPARATION()
    // Preparation has no tasks: this value is ready before route selection.
    route_settings = prepared.route_settings.val

    // These groups follow the shared `take:` order used by route workflows.
    def ranking_parameters = [
        prepared.rankscore_filter,
        prepared.rankscore_softwares,
        prepared.phen2gene_top_n,
        prepared.gnomad_af_ceiling,
        prepared.minimum_genotype_quality,
        prepared.minimum_allele_depth,
        prepared.rankvar_filter
    ]
    def interpretation_parameters = [
        prepared.inheritance_mode,
        prepared.include_clinvar_report,
        prepared.allow_unphased_comphet
    ]
    // Select one nuclear route; mitochondrial evidence is added where supported.
    if (route_settings.route == 'annotated_ngs') {
        PREANNOTATED_ALL_NGS(
            prepared.input_annotated_ngs, prepared.clinical_metadata_by_sample,
            prepared.reference_bundle, prepared.eh_variant_catalog,
            ranking_parameters[0], ranking_parameters[1], ranking_parameters[2], ranking_parameters[3], ranking_parameters[4], ranking_parameters[5], ranking_parameters[6], interpretation_parameters[0], interpretation_parameters[1], interpretation_parameters[2],
            prepared.mitochondrial_reference_bundle, prepared.mito_contig, route_settings.mitochondrial_enabled
        )
    }
    else if (route_settings.route == 'annotated_snv_sv') {
        PREANNOTATED_SV(
            prepared.input_annotated_snv_sv, prepared.clinical_metadata_by_sample,
            ranking_parameters[0], ranking_parameters[1], ranking_parameters[2], ranking_parameters[3], ranking_parameters[4], ranking_parameters[5], ranking_parameters[6], interpretation_parameters[0], interpretation_parameters[1], interpretation_parameters[2]
        )
    }
    else if (route_settings.route == 'annotated_called_ngs') {
        PREANNOTATED_SNV_CALLED_SV_NGS(
            prepared.input_annotated_called_ngs, prepared.clinical_metadata_by_sample,
            prepared.reference_bundle, prepared.reference_bundle, prepared.eh_variant_catalog,
            ranking_parameters[0], ranking_parameters[1], ranking_parameters[2], ranking_parameters[3], ranking_parameters[4], ranking_parameters[5], ranking_parameters[6], interpretation_parameters[0], interpretation_parameters[1], interpretation_parameters[2],
            prepared.mitochondrial_reference_bundle, prepared.mito_contig, route_settings.mitochondrial_enabled
        )
    }
    else if (route_settings.route == 'annotated_snv') {
        PREANNOTATED_SMALL_VARIANT(
            prepared.input_annotated_snv, prepared.clinical_metadata_by_sample,
            ranking_parameters[0], ranking_parameters[1], ranking_parameters[2], ranking_parameters[3], ranking_parameters[4], ranking_parameters[5], ranking_parameters[6], interpretation_parameters[0], interpretation_parameters[1], interpretation_parameters[2]
        )
    }
    else if (route_settings.route == 'vcf') {
        ALIGNMENT_VCF(
            prepared.input_vcf, prepared.clinical_metadata_by_sample,
            prepared.reference_bundle, route_settings.effective_mode,
            ranking_parameters[0], ranking_parameters[1], ranking_parameters[2], ranking_parameters[3], ranking_parameters[4], ranking_parameters[5], ranking_parameters[6],
            route_settings.phenotype_is_clinical_note, route_settings.phenotype_targeting_enabled,
            interpretation_parameters[0], interpretation_parameters[1], interpretation_parameters[2]
        )
    }
    else if (route_settings.route == 'alignment') {
        if (route_settings.sequencing_type == 'short') {
            mito_report_tsv = Channel.empty()
            if (route_settings.mitochondrial_enabled == 'yes') {
                // Mitochondrial callers need only the sample-keyed alignment; phenotype stays nuclear-side.
                short_mito_input = prepared.mito_input_bam.map {
                    sample_id, bam_file, bai_file, phenotype_file -> tuple(sample_id, bam_file, bai_file)
                }
                ALIGNMENT_NGS_MITO(
                    short_mito_input,
                    prepared.mitochondrial_reference_bundle,
                    prepared.mito_contig
                )
                mito_report_tsv = ALIGNMENT_NGS_MITO.out.prioritized_tsv
            }
            ngs_analysis_mode = route_settings.effective_mode ?: 'combined'
            ALIGNMENT_NGS(
                prepared.input_bam, prepared.clinical_metadata_by_sample,
                prepared.reference_bundle, prepared.gatk_reference_bundle, prepared.eh_variant_catalog,
                ngs_analysis_mode,
                ranking_parameters[0], ranking_parameters[1], ranking_parameters[2], ranking_parameters[3], ranking_parameters[4], ranking_parameters[5], ranking_parameters[6],
                route_settings.phenotype_is_clinical_note, route_settings.phenotype_targeting_enabled, route_settings.short_read_small_variant_caller,
                interpretation_parameters[0], interpretation_parameters[1], interpretation_parameters[2],
                mito_report_tsv, route_settings.mitochondrial_enabled
            )
        }
        else { // Long reads
            mito_report_tsv = Channel.empty()
            if (route_settings.mitochondrial_enabled == 'yes') {
                // Mitochondrial callers need only the sample-keyed alignment; phenotype stays nuclear-side.
                long_mito_input = prepared.mito_input_bam.map {
                    sample_id, bam_file, bai_file, phenotype_file -> tuple(sample_id, bam_file, bai_file)
                }
                ALIGNMENT_LONG_MITO(
                    long_mito_input,
                    prepared.reference_bundle,
                    prepared.mito_contig
                )
                mito_report_tsv = ALIGNMENT_LONG_MITO.out.prioritized_tsv
            }
            long_analysis_mode = route_settings.effective_mode ?: 'combined'
            ALIGNMENT_LONG(
                prepared.input_bam, prepared.clinical_metadata_by_sample, route_settings.batch_input,
                prepared.reference_bundle, long_analysis_mode,
                ranking_parameters[0], ranking_parameters[1], ranking_parameters[2], ranking_parameters[3], ranking_parameters[4], ranking_parameters[5], ranking_parameters[6],
                route_settings.phenotype_is_clinical_note, route_settings.phenotype_targeting_enabled, route_settings.long_read_small_variant_caller,
                interpretation_parameters[0], interpretation_parameters[1], interpretation_parameters[2],
                mito_report_tsv, route_settings.mitochondrial_enabled
            )
        }
    }
}
