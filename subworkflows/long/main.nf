include { LONGPHASE_PROCESSING } from '../longphase_processing'
include { annovar } from '../../modules/annovar/'
include { annovar_sv } from '../../modules/annovar_sv/'
include { clair3 } from '../../modules/clair3/'
include { common_sv_filter } from '../../modules/common_sv_filter/'
include { MATERIALIZE_PUBLIC_ARTIFACT } from '../../modules/local/legacy_artifact_compat/'
include { nanorepeat } from '../../modules/nanorepeat/'
include { nanocaller } from '../../modules/nanocaller/'
include { SNIFFLES } from '../../modules/nf-core/sniffles/'
include { phen2gene } from '../../modules/phen2gene/'
include { phenogpt2 } from '../../modules/phenogpt2/'
include { phenosv } from '../../modules/phenosv/'
include { phenotagger } from '../../modules/phenotagger/'
include { rankscore } from '../../modules/rankscore_analysis/'
include { rankvar } from '../../modules/rankvar/'
include { rankvar as rankvar_nanocaller } from '../../modules/rankvar/'
include { reduce_region_phen2gene } from '../../modules/reduce_region_phen2gene/'
include { small_variant_prio } from '../../modules/snp_prio/'
include { survivor } from '../../modules/survivor/'
include { sv_prio } from '../../modules/sv_prio/'
include { variant_html_report_with_repeat; variant_html_report_with_repeat_and_mito } from '../../modules/variant_html_report/'

// Unified long-read dispatcher. Shared sample projections and phenotype
// extraction run once; the SNV, SV, and combined terminal contracts remain
// explicitly guarded so no process receives a placeholder path.
workflow ALIGNMENT_LONG {
    take:
    input_bam
    clinical_metadata_by_sample
    batch_input
    reference_bundle
    analysis_mode
    rankscore_filter
    rankscore_softwares
    phen2gene_top_n
    gnomad_af_ceiling
    minimum_genotype_quality
    minimum_allele_depth
    rankvar_filter
    phenotype_is_clinical_note
    phenotype_targeting_enabled
    small_variant_caller
    inheritance_mode
    include_clinvar_report
    allow_unphased_comphet
    mito_tsv
    mito_mode

    main:
    if (!(analysis_mode in ['snp', 'sv', 'combined'])) {
        error "Unsupported long-read analysis mode '${analysis_mode}'; expected snp, sv, or combined"
    }

    input_bam.multiMap { sample_id, bam_file, bai_file, phenotype_file ->
        phenotype: tuple(sample_id, phenotype_file)
        alignment: tuple(sample_id, bam_file, bai_file)
    }.set { route_inputs }
    phenotype_inputs_by_sample = route_inputs.phenotype
    alignments_by_sample = route_inputs.alignment

    hpo_files_by_sample = phenotype_inputs_by_sample
    if (phenotype_is_clinical_note == 'yes') {
        if (params.phenotype_extractor.toString().trim().toLowerCase() == 'phenogpt2') {
            hpo_files_by_sample = phenogpt2(phenotype_inputs_by_sample)
        }
        else {
            hpo_files_by_sample = phenotagger(phenotype_inputs_by_sample)
        }
    }
    phenotype_metadata = hpo_files_by_sample
        .join(clinical_metadata_by_sample, failOnMismatch: true, failOnDuplicate: true)
        .map { sample_id, hpo_path, age_of_onset -> tuple(sample_id, hpo_path, age_of_onset) }

    if (analysis_mode in ['snp', 'combined']) {
        phen2gene_result = phen2gene(hpo_files_by_sample)
        caller_regions = alignments_by_sample.map { sample_id, bam_file, bai_file -> tuple(sample_id, []) }
        if (phenotype_targeting_enabled == 'yes') {
            phenotype_target_regions = reduce_region_phen2gene(phen2gene_result, reference_bundle, phen2gene_top_n)
            caller_regions = phenotype_target_regions
        }
        caller_input = alignments_by_sample.join(caller_regions, failOnMismatch: true, failOnDuplicate: true)
        if (small_variant_caller == 'nanocaller') {
            small_variant_calls = nanocaller(caller_input, reference_bundle)
        }
        else {
            small_variant_calls = clair3(caller_input, reference_bundle)
        }

        if (phenotype_targeting_enabled == 'yes') {
            annovar_input = small_variant_calls
                .join(phenotype_target_regions, failOnMismatch: true, failOnDuplicate: true)
                .map { sample_id, vcf_file, bed_file -> tuple(sample_id, vcf_file, bed_file) }
        }
        else {
            annovar_input = small_variant_calls.map { sample_id, vcf_file -> tuple(sample_id, vcf_file, []) }
        }
        annovar_for_downstream = annovar(annovar_input)
        annovar_txt_by_sample = annovar_for_downstream.map { sample_id, annovar_txt, annovar_vcf -> tuple(sample_id, annovar_txt) }
        annovar_with_ranked_genes = annovar_txt_by_sample.join(
            phen2gene_result,
            failOnMismatch: true,
            failOnDuplicate: true
        )
        annovar_with_ranked_genes_and_hpo = annovar_with_ranked_genes.join(
            hpo_files_by_sample,
            failOnMismatch: true,
            failOnDuplicate: true
        )
        rankscore_result = rankscore(
            annovar_with_ranked_genes,
            gnomad_af_ceiling,
            rankscore_filter,
            rankscore_softwares,
            minimum_genotype_quality,
            minimum_allele_depth,
            phen2gene_top_n
        )
        if (small_variant_caller == 'nanocaller') {
            rankvar_result = rankvar_nanocaller(
                annovar_with_ranked_genes_and_hpo,
                gnomad_af_ceiling,
                minimum_genotype_quality,
                minimum_allele_depth,
                rankvar_filter,
                'nanocaller',
                params.nanocaller_dp
            )
        }
        else {
            rankvar_result = rankvar(
                annovar_with_ranked_genes_and_hpo,
                gnomad_af_ceiling,
                minimum_genotype_quality,
                minimum_allele_depth,
                rankvar_filter,
                'standard',
                params.nanocaller_dp
            )
        }

        if (analysis_mode == 'snp') {
            snp_annovar_vcf_by_sample = annovar_for_downstream
                .map { sample_id, annovar_txt, annovar_vcf -> tuple(sample_id, annovar_vcf) }
            small_variant_priority_input = rankscore_result
                .join(rankvar_result, failOnMismatch: true, failOnDuplicate: true)
                .join(snp_annovar_vcf_by_sample, failOnMismatch: true, failOnDuplicate: true)
                .join(phenotype_metadata, failOnMismatch: true, failOnDuplicate: true)
            small_variant_prio(
                small_variant_priority_input,
                inheritance_mode,
                include_clinvar_report,
                allow_unphased_comphet
            )
        }
    }

    if (analysis_mode in ['sv', 'combined']) {
        sniffles_input = alignments_by_sample.map { sample_id, bam_file, index_file ->
            tuple([id: sample_id, single_end: false], bam_file, index_file)
        }
        sniffles_ref_fasta = reference_bundle.map { fasta, fai -> tuple([id: 'reference'], fasta) }
        sniffles_tandem_file = reference_bundle.map { fasta, fai -> tuple([id: 'reference'], []) }
        SNIFFLES(sniffles_input, sniffles_ref_fasta, sniffles_tandem_file, true, false)
        sniffles_public_artifact_input = SNIFFLES.out.vcf.map { meta, vcf_file ->
            tuple(meta, vcf_file, 'sniffles', "${meta.id}.sniffles.vcf")
        }
        sniffles_result = MATERIALIZE_PUBLIC_ARTIFACT(sniffles_public_artifact_input)
            .map { meta, vcf_file, source_tool -> tuple(meta.id, vcf_file) }

        sniffles_result_annovar = sniffles_result.map { sample_id, vcf_file -> tuple(sample_id, vcf_file, 'called') }
        annovar_sv_for_downstream = annovar_sv(sniffles_result_annovar)
        if (params.common_sv_filter.toString().trim().toLowerCase() == 'yes') {
            common_sv_filter(annovar_sv_for_downstream)
            annovar_sv_for_downstream = common_sv_filter.out.filtered_vcf
        }
        survivor_result = survivor(annovar_sv_for_downstream)
        phenosv_input = survivor_result.join(hpo_files_by_sample, failOnMismatch: true, failOnDuplicate: true)
        phenosv_result = phenosv(phenosv_input)
        nanorepeat(alignments_by_sample, reference_bundle)

        if (analysis_mode == 'sv') {
            sv_priority_input = phenosv_result
                .join(annovar_sv_for_downstream, failOnMismatch: true, failOnDuplicate: true)
                .join(phenotype_metadata, failOnMismatch: true, failOnDuplicate: true)
            sv_prio(
                sv_priority_input,
                inheritance_mode,
                include_clinvar_report,
                allow_unphased_comphet
            )
        }
    }

    if (analysis_mode == 'combined') {
        longphase_annovar_vcf_by_sample = annovar_for_downstream
            .map { sample_id, annovar_txt, annovar_vcf -> tuple(sample_id, annovar_vcf) }
        join_vcf_bam = longphase_annovar_vcf_by_sample.join(
            alignments_by_sample,
            failOnMismatch: true,
            failOnDuplicate: true
        )
        join_vcf_bam_sv = sniffles_result
            .join(annovar_sv_for_downstream, failOnMismatch: true, failOnDuplicate: true)
            .join(join_vcf_bam, failOnMismatch: true, failOnDuplicate: true)
        join_vcf_bam_phenosv = phenosv_result.join(
            join_vcf_bam_sv,
            failOnMismatch: true,
            failOnDuplicate: true
        )
        join_vcf_bam_rankscore = rankscore_result.join(
            join_vcf_bam_phenosv,
            failOnMismatch: true,
            failOnDuplicate: true
        )
        join_vcf_bam_rankvar = rankvar_result.join(
            join_vcf_bam_rankscore,
            failOnMismatch: true,
            failOnDuplicate: true
        )
        longphase_input = join_vcf_bam_rankvar.join(
            phenotype_metadata,
            failOnMismatch: true,
            failOnDuplicate: true
        )
        longphase_result = LONGPHASE_PROCESSING(
            longphase_input,
            reference_bundle,
            batch_input,
            inheritance_mode,
            include_clinvar_report,
            allow_unphased_comphet
        )

        prio_report_input = longphase_result.prio_vcf
            .join(longphase_result.prio_gene, failOnMismatch: true, failOnDuplicate: true)
            .join(nanorepeat.out, failOnMismatch: true, failOnDuplicate: true)
            .map { sample_id, prio_vcf, prio_gene_report, repeat_tsv ->
                tuple(sample_id, prio_vcf, prio_gene_report, repeat_tsv)
            }
        if (mito_mode == 'yes') {
            prio_report_input_with_mito = prio_report_input
                .join(
                    mito_tsv.map { sample_id, mito_report -> tuple(sample_id, mito_report) },
                    failOnMismatch: true,
                    failOnDuplicate: true
                )
                .map { sample_id, prio_vcf, prio_gene_report, repeat_tsv, mito_report ->
                    tuple(sample_id, prio_vcf, prio_gene_report, repeat_tsv, mito_report)
                }
            variant_html_report_with_repeat_and_mito(prio_report_input_with_mito)
        }
        else {
            variant_html_report_with_repeat(prio_report_input)
        }
    }
}
