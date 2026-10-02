include { PREANNOTATED_SNV_EVIDENCE } from '../preannotated_shared/main'
include { small_variant_prio } from '../../modules/snp_prio/'

// Shared core: pre-annotated ANNOVAR TXT/VCF SNP prioritization.
workflow PREANNOTATED_SMALL_VARIANT {
    take:
    input_annotated_snv
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
    PREANNOTATED_SNV_EVIDENCE(
        input_annotated_snv,
        clinical_metadata_by_sample,
        rankscore_filter,
        rankscore_softwares,
        phen2gene_top_n,
        gnomad_af_ceiling,
        minimum_genotype_quality,
        minimum_allele_depth,
        rankvar_filter
    )
    rankscore_with_rankvar = PREANNOTATED_SNV_EVIDENCE.out.rankscore.join(
        PREANNOTATED_SNV_EVIDENCE.out.rankvar,
        failOnMismatch: true,
        failOnDuplicate: true
    )
    annovar_vcf_for_prio = PREANNOTATED_SNV_EVIDENCE.out.annovar.map {
        sample_id, annovar_txt, annovar_vcf -> tuple(sample_id, annovar_vcf)
    }
    small_variant_priority_input = rankscore_with_rankvar.join(
        annovar_vcf_for_prio,
        failOnMismatch: true,
        failOnDuplicate: true
    )
    small_variant_priority_with_phenotype = small_variant_priority_input.join(
        PREANNOTATED_SNV_EVIDENCE.out.phenotype_metadata,
        failOnMismatch: true,
        failOnDuplicate: true
    )
    small_variant_prio(small_variant_priority_with_phenotype, inheritance_mode, include_clinvar_report, allow_unphased_comphet)
    prio_vcf_result = small_variant_prio.out.prio_vcf
    prio_gene_vcf_result = small_variant_prio.out.prio_gene_vcf

    emit:
    validated_annovar = PREANNOTATED_SNV_EVIDENCE.out.annovar
    phen2gene = PREANNOTATED_SNV_EVIDENCE.out.phen2gene
    rankscore = PREANNOTATED_SNV_EVIDENCE.out.rankscore
    rankvar = PREANNOTATED_SNV_EVIDENCE.out.rankvar
    prio_vcf = prio_vcf_result
    prio_gene_vcf = prio_gene_vcf_result
}
