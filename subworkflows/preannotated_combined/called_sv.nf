include { PREANNOTATED_SNV_EVIDENCE; PREANNOTATED_COMBINED_EVIDENCE } from '../preannotated_shared/main'
include { EXPANSIONHUNTER } from '../../modules/nf-core/expansionhunter/'
include { eh_filter } from '../../modules/eh_filter/'
include { MATERIALIZE_PUBLIC_ARTIFACT as MANTA_PUBLIC_ARTIFACT } from '../../modules/local/legacy_artifact_compat/'
include { MATERIALIZE_PUBLIC_ARTIFACT as EXPANSIONHUNTER_PUBLIC_ARTIFACT } from '../../modules/local/legacy_artifact_compat/'
include { ngs_prio } from '../../modules/ngs_prio/'
include { annovar_sv } from '../../modules/annovar_sv/'
include { common_sv_filter } from '../../modules/common_sv_filter/'
include { survivor } from '../../modules/survivor/'
include { phenosv } from '../../modules/phenosv/'
include { MANTA_GERMLINE } from '../../modules/nf-core/manta/germline/'
include { xtea } from '../../modules/xtea/'
include { normalize_shortread_alignment } from '../../modules/normalize_shortread_alignment/'
include { cnvnator } from '../../modules/cnvnator/'
include { truvari_shortread_sv_merge } from '../../modules/truvari_shortread_sv_merge/'
include { ALIGNMENT_NGS_MITO } from '../mitochondrial/ngs'
include { variant_html_report_with_repeat; variant_html_report_with_repeat_and_mito } from '../../modules/variant_html_report/'

// Sample-keyed route: imported ANNOVAR SNV + called short-read SV/CNV/(optional) mito analysis.
workflow PREANNOTATED_SNV_CALLED_SV_NGS {
    take:
    input_annotated_called_ngs
    clinical_metadata_by_sample
    reference_bundle
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
    canonical_snv_input = input_annotated_called_ngs.map {
        sample_id, annovar_txt, annovar_vcf, bam_file, bai_file, phenotype_path, phenotype_format ->
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

    alignments_by_sample = input_annotated_called_ngs.map { sample_id, annovar_txt, annovar_vcf, bam_file, bai_file, phenotype_path, phenotype_format ->
        tuple(sample_id, bam_file, bai_file)
    }
    all_alignments_by_sample = input_annotated_called_ngs.map { sample_id, annovar_txt, annovar_vcf, bam_file, bai_file, phenotype_path, phenotype_format ->
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
    expansionhunter_json = EXPANSIONHUNTER_PUBLIC_ARTIFACT(expansionhunter_public_artifact_input).map { meta, json_file, source_tool -> tuple(meta.id, json_file) }
    eh_filter(expansionhunter_json)

    manta_input = all_alignments_by_sample.map { sample_id, bam_file, index_file -> tuple([id: sample_id, single_end: false], bam_file, index_file, [], []) }
    manta_ref_fasta = reference_bundle.map { fasta, fai -> tuple([id: 'reference'], fasta) }
    manta_ref_fai = reference_bundle.map { fasta, fai -> tuple([id: 'reference'], fai) }
    manta_config = reference_bundle.map { fasta, fai -> [] }
    MANTA_GERMLINE(manta_input, manta_ref_fasta, manta_ref_fai, manta_config)
    manta_public_artifact_input = MANTA_GERMLINE.out.diploid_sv_vcf.map { meta, vcf_file -> tuple(meta, vcf_file, 'manta', "${meta.id}_manta.vcf") }
    manta_result = MANTA_PUBLIC_ARTIFACT(manta_public_artifact_input).map { meta, vcf_file, source_tool -> tuple(meta.id, vcf_file) }
    xtea_mode = params.xtea ? params.xtea.toString().trim().toLowerCase() : "no"
    xtea_vcf = Channel.empty()
    if ( xtea_mode == "yes" ) {
        xtea_input = all_alignments_by_sample.map { sample_id, bam_file, index_file ->
            tuple([id: sample_id], bam_file, index_file)
        }
        xtea(xtea_input, reference_bundle)
        xtea_vcf = xtea.out.vcf.map { meta, vcf -> tuple(meta.id, vcf) }
    }

    cnvnator_mode = params.cnvnator ? params.cnvnator.toString().trim().toLowerCase() : "yes"
    // Optional callers produce separate per-sample channels only when enabled.
    // Each branch therefore joins exactly the available callers before Truvari.
    if ( cnvnator_mode != "no" && xtea_mode == "yes" ) {
        normalized_bam = normalize_shortread_alignment(all_alignments_by_sample, reference_bundle)
        cnvnator(normalized_bam, reference_bundle, params.cnvnator_bin_size)
        merged_sv_input = manta_result.join(cnvnator.out.vcf, failOnMismatch: true, failOnDuplicate: true).join(xtea_vcf, failOnMismatch: true, failOnDuplicate: true).map { sample_id, manta_vcf, cnvnator_vcf, xtea_vcf ->
            tuple(sample_id, [manta_vcf, cnvnator_vcf, xtea_vcf])
        }
        truvari_shortread_sv_merge(merged_sv_input, reference_bundle)
        sv_result = truvari_shortread_sv_merge.out.merged_vcf
    }
    else if ( cnvnator_mode != "no" ) {
        normalized_bam = normalize_shortread_alignment(all_alignments_by_sample, reference_bundle)
        cnvnator(normalized_bam, reference_bundle, params.cnvnator_bin_size)
        merged_sv_input = manta_result.join(cnvnator.out.vcf, failOnMismatch: true, failOnDuplicate: true).map { sample_id, manta_vcf, cnvnator_vcf ->
            tuple(sample_id, [manta_vcf, cnvnator_vcf])
        }
        truvari_shortread_sv_merge(merged_sv_input, reference_bundle)
        sv_result = truvari_shortread_sv_merge.out.merged_vcf
    }
    else if ( xtea_mode == "yes" ) {
        merged_sv_input = manta_result.join(xtea_vcf, failOnMismatch: true, failOnDuplicate: true).map { sample_id, manta_vcf, xtea_vcf ->
            tuple(sample_id, [manta_vcf, xtea_vcf])
        }
        truvari_shortread_sv_merge(merged_sv_input, reference_bundle)
        sv_result = truvari_shortread_sv_merge.out.merged_vcf
    }
    else {
        sv_result = manta_result
    }

    sv_for_annotation = sv_result
    sv_result_annovar = sv_for_annotation.map { sample_id, vcf_file ->
        tuple(sample_id, vcf_file, "called")
    }
    annovar_sv_result = annovar_sv(sv_result_annovar)
    annovar_sv_for_downstream = annovar_sv_result
    if ( params.common_sv_filter.toString().trim().toLowerCase() == "yes" ) {
        common_sv_filter(annovar_sv_for_downstream)
        annovar_sv_for_downstream = common_sv_filter.out.filtered_vcf
    }
    survivor_result = survivor(annovar_sv_for_downstream)
    phenosv_input = survivor_result.join(hpo_paths, failOnMismatch: true, failOnDuplicate: true)
    phenosv_result = phenosv(phenosv_input)

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
