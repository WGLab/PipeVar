include { PREANNOTATED_SNV_EVIDENCE; PREANNOTATED_IMPORTED_SV_CHAIN; PREANNOTATED_COMBINED_EVIDENCE } from '../preannotated_shared/main'
include { EXPANSIONHUNTER } from '../../modules/nf-core/expansionhunter/'
include { eh_filter } from '../../modules/eh_filter/'
include { MATERIALIZE_PUBLIC_ARTIFACT } from '../../modules/local/legacy_artifact_compat/'
include { ngs_prio } from '../../modules/ngs_prio/'
include { ALIGNMENT_NGS_MITO } from '../mitochondrial/ngs'
include { variant_html_report_with_repeat; variant_html_report_with_repeat_and_mito } from '../../modules/variant_html_report/'

// Sample-keyed route: imported ANNOVAR SNV + short-read SV/CNV/(optional) mito analysis.
workflow PREANNOTATED_ALL_NGS {
    take:
    input_annotated_ngs
    clinical_metadata_by_sample
    eh_ref_fa
    eh_variant_catalog
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
    mitochondrial_reference_bundle
    mito_contig
    mito_mode

    main:
    canonical_snv_input = input_annotated_ngs.map {
        sample_id, annovar_txt, annovar_vcf, annovar_sv_vcf, bam_file, bai_file, phenotype_path, phenotype_format ->
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
    annovar_for_downstream = PREANNOTATED_SNV_EVIDENCE.out.annovar
    hpo_paths = PREANNOTATED_SNV_EVIDENCE.out.hpo
    rankscore_result = PREANNOTATED_SNV_EVIDENCE.out.rankscore
    rankvar_result = PREANNOTATED_SNV_EVIDENCE.out.rankvar
    phenotype_metadata = PREANNOTATED_SNV_EVIDENCE.out.phenotype_metadata

    alignments_by_sample = input_annotated_ngs.map { sample_id, annovar_txt, annovar_vcf, annovar_sv_vcf, bam_file, bai_file, phenotype_path, phenotype_format ->
        tuple(sample_id, bam_file, bai_file)
    }
    expansionhunter_input = alignments_by_sample.map { sample_id, bam_file, bai_file -> tuple([id: sample_id, single_end: false], bam_file, bai_file) }
    expansionhunter_ref_fasta = eh_ref_fa.map { fasta, fai -> tuple([id: 'reference'], fasta) }
    expansionhunter_ref_fai = eh_ref_fa.map { fasta, fai -> tuple([id: 'reference'], fai) }
    expansionhunter_catalog = eh_variant_catalog.map { catalog -> tuple([id: 'expansionhunter_catalog'], catalog) }
    EXPANSIONHUNTER(
        expansionhunter_input,
        expansionhunter_ref_fasta,
        expansionhunter_ref_fai,
        expansionhunter_catalog
    )
    expansionhunter_public_artifact_input = EXPANSIONHUNTER.out.json.map { meta, json_file -> tuple(meta, json_file, 'expansionhunter', "${meta.id}.json") }
    expansionhunter_json = MATERIALIZE_PUBLIC_ARTIFACT(expansionhunter_public_artifact_input).map { meta, json_file, source_tool -> tuple(meta.id, json_file) }
    eh_filter(expansionhunter_json)

    annotated_sv_input = input_annotated_ngs.map { sample_id, annovar_txt, annovar_vcf, annovar_sv_vcf, bam_file, bai_file, phenotype_path, phenotype_format ->
        tuple(sample_id, annovar_sv_vcf)
    }
    PREANNOTATED_IMPORTED_SV_CHAIN(annotated_sv_input, hpo_paths)
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

    if ( mito_mode == "yes" ) {
        ALIGNMENT_NGS_MITO(alignments_by_sample, mitochondrial_reference_bundle, mito_contig)
    }

    prio_report_input = prio_vcf_result
        .join(prio_gene_result, failOnMismatch: true, failOnDuplicate: true)
        .join(eh_filter.out, failOnMismatch: true, failOnDuplicate: true)
        .map { sample_id, prio_vcf, prio_gene_report, repeat_tsv ->
            tuple(sample_id, prio_vcf, prio_gene_report, repeat_tsv)
        }

    if ( mito_mode == "yes" ) {
        prio_report_input_with_mito = prio_report_input
            .join(ALIGNMENT_NGS_MITO.out.prioritized_tsv, failOnMismatch: true, failOnDuplicate: true)
            .map { sample_id, prio_vcf, prio_gene_report, repeat_tsv, mito_tsv ->
                tuple(sample_id, prio_vcf, prio_gene_report, repeat_tsv, mito_tsv)
            }
        variant_html_report_with_repeat_and_mito(prio_report_input_with_mito)
    }
    else {
        variant_html_report_with_repeat(prio_report_input)
    }

    emit:
    prio_vcf = prio_vcf_result.map { sample_id, vcf -> vcf }
    prio_gene_vcf = prio_gene_result.map { sample_id, vcf -> vcf }
    mito_vcf = mito_mode == "yes" ? ALIGNMENT_NGS_MITO.out.mito_vcf : Channel.empty()
    mito_annotated_tsv = mito_mode == "yes" ? ALIGNMENT_NGS_MITO.out.annotated_tsv : Channel.empty()
    mito_prioritized_tsv = mito_mode == "yes" ? ALIGNMENT_NGS_MITO.out.prioritized_tsv : Channel.empty()
}
