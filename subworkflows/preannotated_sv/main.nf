include { PREANNOTATED_SNV_EVIDENCE; PREANNOTATED_IMPORTED_SV_CHAIN; PREANNOTATED_COMBINED_EVIDENCE } from '../preannotated_shared/main'
include { ngs_prio } from '../../modules/ngs_prio/'
include { variant_html_report_no_repeat } from '../../modules/variant_html_report/'

// Sample-keyed route: imported ANNOVAR SNV and SV inputs without BAM/CRAM-dependent analysis.
workflow PREANNOTATED_SV {
    take:
    input_annotated_snv_sv
    clinical_metadata_by_sample
    rankscore_filter
    rankscore_softwares
    phen2gene_top_n
    gnomad_af_ceiling
    minimum_genotype_quality
    minimum_allele_depth
    rankvar_filter
    inheritance_mode
    include_clinvar_report
    allow_unphased_comphet

    main:
    canonical_snv_input = input_annotated_snv_sv.map {
        sample_id, annovar_txt, annovar_vcf, annovar_sv_vcf, phenotype_path, phenotype_format ->
            tuple(sample_id, annovar_txt, annovar_vcf, phenotype_path, phenotype_format)
    }
    PREANNOTATED_SNV_EVIDENCE(
        canonical_snv_input,
        clinical_metadata_by_sample,
        rankscore_filter,
        rankscore_softwares,
        phen2gene_top_n,
        gnomad_af_ceiling,
        minimum_genotype_quality,
        minimum_allele_depth,
        rankvar_filter
    )
    annotated_sv_input = input_annotated_snv_sv.map { sample_id, annovar_txt, annovar_vcf, annovar_sv_vcf, phenotype_path, phenotype_format ->
        tuple(sample_id, annovar_sv_vcf)
    }
    PREANNOTATED_IMPORTED_SV_CHAIN(annotated_sv_input, PREANNOTATED_SNV_EVIDENCE.out.hpo)
    annovar_for_downstream = PREANNOTATED_SNV_EVIDENCE.out.annovar
    rankscore_result = PREANNOTATED_SNV_EVIDENCE.out.rankscore
    rankvar_result = PREANNOTATED_SNV_EVIDENCE.out.rankvar
    phenotype_metadata = PREANNOTATED_SNV_EVIDENCE.out.phenotype_metadata
    annovar_sv_for_downstream = PREANNOTATED_IMPORTED_SV_CHAIN.out.annovar_sv
    phenosv_result = PREANNOTATED_IMPORTED_SV_CHAIN.out.phenosv

    PREANNOTATED_COMBINED_EVIDENCE(
        annovar_for_downstream,
        rankscore_result,
        rankvar_result,
        phenosv_result,
        annovar_sv_for_downstream,
        phenotype_metadata
    )
    ngs_prio(PREANNOTATED_COMBINED_EVIDENCE.out.priority_input, inheritance_mode, include_clinvar_report, allow_unphased_comphet)
    prio_vcf_result = ngs_prio.out.prio_vcf
    prio_gene_result = ngs_prio.out.prio_gene_vcf

    prio_report_input = prio_vcf_result
        .join(prio_gene_result, failOnMismatch: true, failOnDuplicate: true)
    variant_html_report_no_repeat(prio_report_input)

    emit:
    prio_vcf = prio_vcf_result.map { sample_id, vcf -> vcf }
    prio_gene_vcf = prio_gene_result.map { sample_id, vcf -> vcf }
    html_report = variant_html_report_no_repeat.out
}
