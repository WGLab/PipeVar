include { annovar } from '../../modules/annovar/'
include { annovar_sv } from '../../modules/annovar_sv/'
include { common_sv_filter } from '../../modules/common_sv_filter/'
include { phen2gene } from '../../modules/phen2gene/'
include { phenogpt2 } from '../../modules/phenogpt2/'
include { phenosv } from '../../modules/phenosv/'
include { phenotagger } from '../../modules/phenotagger/'
include { rankscore } from '../../modules/rankscore_analysis/'
include { rankvar } from '../../modules/rankvar/'
include { reduce_region_phen2gene } from '../../modules/reduce_region_phen2gene/'
include { small_variant_prio } from '../../modules/snp_prio/'
include { survivor } from '../../modules/survivor/'
include { sv_prio } from '../../modules/sv_prio/'

// Sample-keyed raw-VCF route. The public CLI requires either SNP or SV mode;
// phenotype preparation is shared while annotation and prioritization remain
// mode-specific because their biological and file contracts differ.
workflow ALIGNMENT_VCF {
    take:
    input_vcf
    clinical_metadata_by_sample
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
    inheritance_mode
    include_clinvar_report
    allow_unphased_comphet

    main:
    if (!(analysis_mode in ['snp', 'sv'])) {
        error "Unsupported VCF analysis mode '${analysis_mode}'; expected snp or sv"
    }

    phenotype_inputs_by_sample = input_vcf.map { sample_id, vcf_file, phenotype_file ->
        tuple(sample_id, phenotype_file)
    }
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

    if (analysis_mode == 'snp') {
        snv_for_annotation = input_vcf.map { sample_id, vcf_file, phenotype_file -> tuple(sample_id, vcf_file) }
        phen2gene_result = phen2gene(hpo_files_by_sample)

        if (phenotype_targeting_enabled == 'yes') {
            phenotype_target_regions = reduce_region_phen2gene(phen2gene_result, reference_bundle, phen2gene_top_n)
            annovar_input = snv_for_annotation
                .join(phenotype_target_regions, failOnMismatch: true, failOnDuplicate: true)
                .map { sample_id, vcf_file, bed_file -> tuple(sample_id, vcf_file, bed_file) }
        }
        else {
            annovar_input = snv_for_annotation.map { sample_id, vcf_file -> tuple(sample_id, vcf_file, []) }
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
        rankvar_result = rankvar(
            annovar_with_ranked_genes_and_hpo,
            gnomad_af_ceiling,
            minimum_genotype_quality,
            minimum_allele_depth,
            rankvar_filter,
            'standard',
            params.nanocaller_dp
        )

        rankscore_with_rankvar = rankscore_result.join(rankvar_result, failOnMismatch: true, failOnDuplicate: true)
        annovar_vcf_by_sample = annovar_for_downstream.map { sample_id, annovar_txt, annovar_vcf -> tuple(sample_id, annovar_vcf) }
        small_variant_priority_input = rankscore_with_rankvar
            .join(annovar_vcf_by_sample, failOnMismatch: true, failOnDuplicate: true)
            .join(phenotype_metadata, failOnMismatch: true, failOnDuplicate: true)
        small_variant_prio(
            small_variant_priority_input,
            inheritance_mode,
            include_clinvar_report,
            allow_unphased_comphet
        )
    }

    if (analysis_mode == 'sv') {
        sv_for_annotation = input_vcf.map { sample_id, vcf_file, phenotype_file -> tuple(sample_id, vcf_file, 'called') }
        annovar_sv_for_downstream = annovar_sv(sv_for_annotation)
        if (params.common_sv_filter.toString().trim().toLowerCase() == 'yes') {
            common_sv_filter(annovar_sv_for_downstream)
            annovar_sv_for_downstream = common_sv_filter.out.filtered_vcf
        }
        survivor_result = survivor(annovar_sv_for_downstream)
        phenosv_input = survivor_result.join(hpo_files_by_sample, failOnMismatch: true, failOnDuplicate: true)
        phenosv_result = phenosv(phenosv_input)
        sv_prio_input = phenosv_result
            .join(annovar_sv_for_downstream, failOnMismatch: true, failOnDuplicate: true)
            .join(phenotype_metadata, failOnMismatch: true, failOnDuplicate: true)
        sv_prio(
            sv_prio_input,
            inheritance_mode,
            include_clinvar_report,
            allow_unphased_comphet
        )
    }
}
