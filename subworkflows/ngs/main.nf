include { annovar } from '../../modules/annovar/'
include { phen2gene } from '../../modules/phen2gene/'
include { rankscore } from '../../modules/rankscore_analysis/'
include { annovar_sv } from '../../modules/annovar_sv/'
include { common_sv_filter } from '../../modules/common_sv_filter/'
include { survivor } from '../../modules/survivor/'
include { phenosv } from '../../modules/phenosv/'
include { rankvar } from '../../modules/rankvar/'
include { MANTA_GERMLINE } from '../../modules/nf-core/manta/germline/'
include { xtea } from '../../modules/xtea/'
include { normalize_shortread_alignment } from '../../modules/normalize_shortread_alignment/'
include { cnvnator } from '../../modules/cnvnator/'
include { truvari_shortread_sv_merge } from '../../modules/truvari_shortread_sv_merge/'
include { EXPANSIONHUNTER } from '../../modules/nf-core/expansionhunter/'
include { eh_filter } from '../../modules/eh_filter/'
include { DEEPVARIANT_RUNDEEPVARIANT } from '../../modules/nf-core/deepvariant/rundeepvariant/'
include { MATERIALIZE_PUBLIC_ARTIFACT as MANTA_PUBLIC_ARTIFACT } from '../../modules/local/legacy_artifact_compat/'
include { MATERIALIZE_PUBLIC_ARTIFACT as EXPANSIONHUNTER_PUBLIC_ARTIFACT } from '../../modules/local/legacy_artifact_compat/'
include { GATK_SMALL_VARIANT_CALLING } from '../gatk_snp_calling'
include { phenotagger } from '../../modules/phenotagger/'
include { phenogpt2 } from '../../modules/phenogpt2/'
include { reduce_region_phen2gene } from '../../modules/reduce_region_phen2gene/'
include { small_variant_prio } from '../../modules/snp_prio/'
include { sv_prio } from '../../modules/sv_prio/'
include { ngs_prio } from '../../modules/ngs_prio/'
include { variant_html_report_with_repeat; variant_html_report_with_repeat_and_mito } from '../../modules/variant_html_report/'

// Unified short-read route. analysis_mode is one of: snp, sv, combined.
// Input alignment records are (sample_id, alignment, index, phenotype_file).
workflow ALIGNMENT_NGS {
    take:
    input_bam
    clinical_metadata_by_sample
    reference_bundle
    gatk_reference_bundle
    eh_variant_catalog
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
        error "ERROR: Invalid internal short-read analysis mode '${analysis_mode}'."
    }

    phenotype_inputs_by_sample = input_bam.map { sample_id, bam_file, bai_file, phenotype_file ->
        tuple(sample_id, phenotype_file)
    }
    alignments_by_sample = input_bam.map { sample_id, bam_file, bai_file, phenotype_file ->
        tuple(sample_id, bam_file, bai_file)
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
    phenotype_metadata_by_sample = hpo_files_by_sample
        .join(clinical_metadata_by_sample, failOnMismatch: true, failOnDuplicate: true)
        .map { sample_id, hpo_path, age_of_onset -> tuple(sample_id, hpo_path, age_of_onset) }

    // Repeat expansion screening is part of every short-read alignment route.
    expansionhunter_input = alignments_by_sample.map { sample_id, bam_file, bai_file ->
        tuple([id: sample_id, single_end: false], bam_file, bai_file)
    }
    expansionhunter_ref_fasta = reference_bundle.map { fasta, fai -> tuple([id: 'reference'], fasta) }
    expansionhunter_ref_fai = reference_bundle.map { fasta, fai -> tuple([id: 'reference'], fai) }
    expansionhunter_catalog = eh_variant_catalog.map { catalog -> tuple([id: 'expansionhunter_catalog'], catalog) }
    EXPANSIONHUNTER(
        expansionhunter_input,
        expansionhunter_ref_fasta,
        expansionhunter_ref_fai,
        expansionhunter_catalog
    )
    expansionhunter_public_artifact_input = EXPANSIONHUNTER.out.json.map { meta, json_file ->
        tuple(meta, json_file, 'expansionhunter', "${meta.id}.json")
    }
    expansionhunter_json = EXPANSIONHUNTER_PUBLIC_ARTIFACT(expansionhunter_public_artifact_input).map {
        meta, json_file, source_tool -> tuple(meta.id, json_file)
    }
    eh_filter(expansionhunter_json)

    if (analysis_mode in ['snp', 'combined']) {
        phen2gene_result = phen2gene(hpo_files_by_sample)
        caller_regions = alignments_by_sample.map { sample_id, bam_file, bai_file -> tuple(sample_id, []) }
        if (phenotype_targeting_enabled == 'yes') {
            phenotype_target_regions = reduce_region_phen2gene(
                phen2gene_result,
                reference_bundle,
                phen2gene_top_n
            )
            caller_regions = phenotype_target_regions
        }

        if (small_variant_caller == 'haplotypecaller') {
            haplotypecaller_input = alignments_by_sample.join(
                caller_regions,
                failOnMismatch: true,
                failOnDuplicate: true
            )
            gatk_result = GATK_SMALL_VARIANT_CALLING(haplotypecaller_input, gatk_reference_bundle)
            small_variant_calls = gatk_result.vcf
        }
        else {
            deepvariant_input = alignments_by_sample
                .join(caller_regions, failOnMismatch: true, failOnDuplicate: true)
                .map { sample_id, bam_file, bai_file, regions ->
                    tuple([id: sample_id, single_end: false], bam_file, bai_file, regions)
                }
            deepvariant_ref_fasta = reference_bundle.map { fasta, fai -> tuple([id: 'reference'], fasta) }
            deepvariant_ref_fai = reference_bundle.map { fasta, fai -> tuple([id: 'reference'], fai) }
            deepvariant_ref_gzi = reference_bundle.map { fasta, fai -> tuple([id: 'reference'], []) }
            deepvariant_par_bed = reference_bundle.map { fasta, fai -> tuple([id: 'reference'], []) }
            DEEPVARIANT_RUNDEEPVARIANT(
                deepvariant_input,
                deepvariant_ref_fasta,
                deepvariant_ref_fai,
                deepvariant_ref_gzi,
                deepvariant_par_bed
            )
            small_variant_calls = DEEPVARIANT_RUNDEEPVARIANT.out.vcf.map { meta, vcf_file ->
                tuple(meta.id, vcf_file)
            }
        }

        if (phenotype_targeting_enabled == 'yes') {
            annovar_input = small_variant_calls
                .join(phenotype_target_regions, failOnMismatch: true, failOnDuplicate: true)
                .map { sample_id, vcf_file, bed_file -> tuple(sample_id, vcf_file, bed_file) }
        }
        else {
            annovar_input = small_variant_calls.map { sample_id, vcf_file -> tuple(sample_id, vcf_file, []) }
        }
        annovar_result = annovar(annovar_input)
        annovar_txt_by_sample = annovar_result.map { sample_id, annovar_txt, annovar_vcf ->
            tuple(sample_id, annovar_txt)
        }
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
        annovar_vcf_by_sample = annovar_result.map { sample_id, annovar_txt, annovar_vcf ->
            tuple(sample_id, annovar_vcf)
        }

        if (analysis_mode == 'snp') {
            small_variant_priority_input = rankscore_result
                .join(rankvar_result, failOnMismatch: true, failOnDuplicate: true)
                .join(annovar_vcf_by_sample, failOnMismatch: true, failOnDuplicate: true)
            small_variant_priority_with_phenotype = small_variant_priority_input.join(
                phenotype_metadata_by_sample,
                failOnMismatch: true,
                failOnDuplicate: true
            )
            small_variant_prio(
                small_variant_priority_with_phenotype,
                inheritance_mode,
                include_clinvar_report,
                allow_unphased_comphet
            )
        }
    }

    if (analysis_mode in ['sv', 'combined']) {
        manta_input = alignments_by_sample.map { sample_id, bam_file, index_file ->
            tuple([id: sample_id, single_end: false], bam_file, index_file, [], [])
        }
        manta_ref_fasta = reference_bundle.map { fasta, fai -> tuple([id: 'reference'], fasta) }
        manta_ref_fai = reference_bundle.map { fasta, fai -> tuple([id: 'reference'], fai) }
        manta_config = reference_bundle.map { fasta, fai -> [] }
        MANTA_GERMLINE(manta_input, manta_ref_fasta, manta_ref_fai, manta_config)
        manta_public_artifact_input = MANTA_GERMLINE.out.diploid_sv_vcf.map { meta, vcf_file ->
            tuple(meta, vcf_file, 'manta', "${meta.id}_manta.vcf")
        }
        manta_result = MANTA_PUBLIC_ARTIFACT(manta_public_artifact_input).map {
            meta, vcf_file, source_tool -> tuple(meta.id, vcf_file)
        }

        xtea_mode = params.xtea ? params.xtea.toString().trim().toLowerCase() : 'no'
        xtea_vcf = null
        if (xtea_mode == 'yes') {
            xtea_input = alignments_by_sample.map { sample_id, bam_file, index_file ->
                tuple([id: sample_id], bam_file, index_file)
            }
            xtea(xtea_input, reference_bundle)
            xtea_vcf = xtea.out.vcf.map { meta, vcf_file -> tuple(meta.id, vcf_file) }
        }

        cnvnator_mode = params.cnvnator ? params.cnvnator.toString().trim().toLowerCase() : 'yes'
        if (cnvnator_mode != 'no' && xtea_mode == 'yes') {
            normalized_bam = normalize_shortread_alignment(alignments_by_sample, reference_bundle)
            cnvnator(normalized_bam, reference_bundle, params.cnvnator_bin_size)
            merged_sv_input = manta_result
                .join(cnvnator.out.vcf, failOnMismatch: true, failOnDuplicate: true)
                .join(xtea_vcf, failOnMismatch: true, failOnDuplicate: true)
                .map { sample_id, manta_vcf, cnvnator_vcf, xtea_vcf_file ->
                    tuple(sample_id, [manta_vcf, cnvnator_vcf, xtea_vcf_file])
                }
            truvari_shortread_sv_merge(merged_sv_input, reference_bundle)
            sv_result = truvari_shortread_sv_merge.out.merged_vcf
        }
        else if (cnvnator_mode != 'no') {
            normalized_bam = normalize_shortread_alignment(alignments_by_sample, reference_bundle)
            cnvnator(normalized_bam, reference_bundle, params.cnvnator_bin_size)
            merged_sv_input = manta_result
                .join(cnvnator.out.vcf, failOnMismatch: true, failOnDuplicate: true)
                .map { sample_id, manta_vcf, cnvnator_vcf -> tuple(sample_id, [manta_vcf, cnvnator_vcf]) }
            truvari_shortread_sv_merge(merged_sv_input, reference_bundle)
            sv_result = truvari_shortread_sv_merge.out.merged_vcf
        }
        else if (xtea_mode == 'yes') {
            merged_sv_input = manta_result
                .join(xtea_vcf, failOnMismatch: true, failOnDuplicate: true)
                .map { sample_id, manta_vcf, xtea_vcf_file -> tuple(sample_id, [manta_vcf, xtea_vcf_file]) }
            truvari_shortread_sv_merge(merged_sv_input, reference_bundle)
            sv_result = truvari_shortread_sv_merge.out.merged_vcf
        }
        else {
            sv_result = manta_result
        }

        sv_result_annovar = sv_result.map { sample_id, vcf_file -> tuple(sample_id, vcf_file, 'called') }
        annovar_sv_result = annovar_sv(sv_result_annovar)
        annovar_sv_for_downstream = annovar_sv_result
        if (params.common_sv_filter.toString().trim().toLowerCase() == 'yes') {
            common_sv_filter(annovar_sv_for_downstream)
            annovar_sv_for_downstream = common_sv_filter.out.filtered_vcf
        }
        survivor_result = survivor(annovar_sv_for_downstream)
        phenosv_input = survivor_result.join(
            hpo_files_by_sample,
            failOnMismatch: true,
            failOnDuplicate: true
        )
        phenosv_result = phenosv(phenosv_input)

        if (analysis_mode == 'sv') {
            sv_priority_input = phenosv_result.join(
                annovar_sv_for_downstream,
                failOnMismatch: true,
                failOnDuplicate: true
            )
            sv_priority_with_phenotype = sv_priority_input.join(
                phenotype_metadata_by_sample,
                failOnMismatch: true,
                failOnDuplicate: true
            )
            sv_prio(
                sv_priority_with_phenotype,
                inheritance_mode,
                include_clinvar_report,
                allow_unphased_comphet
            )
        }
    }

    if (analysis_mode == 'combined') {
        phenosv_annovar_snv = phenosv_result.join(
            annovar_vcf_by_sample,
            failOnMismatch: true,
            failOnDuplicate: true
        )
        sv_join = phenosv_annovar_snv.join(
            annovar_sv_for_downstream,
            failOnMismatch: true,
            failOnDuplicate: true
        )
        rankscore_join = sv_join.join(rankscore_result, failOnMismatch: true, failOnDuplicate: true)
        rankvar_join = rankscore_join.join(rankvar_result, failOnMismatch: true, failOnDuplicate: true)
        combined_evidence_by_sample = rankvar_join
            .join(phenotype_metadata_by_sample, failOnMismatch: true, failOnDuplicate: true)
            .map {
                sample_id, sv_phenosv_evidence, snv_vcf, sv_vcf, snv_rankscore,
                snv_clinvar_evidence, snv_rankvar, hpo_path, age_of_onset ->
                tuple(sample_id, [
                    snv_rankvar: snv_rankvar,
                    snv_rankscore: snv_rankscore,
                    snv_clinvar_evidence: snv_clinvar_evidence,
                    sv_phenosv_evidence: sv_phenosv_evidence,
                    sv_vcf: sv_vcf,
                    snv_vcf: snv_vcf,
                    hpo_path: hpo_path,
                    age_of_onset: age_of_onset
                ])
            }
        combined_priority_input = combined_evidence_by_sample.map { sample_id, evidence ->
            tuple(
                sample_id,
                evidence.snv_rankvar,
                evidence.snv_rankscore,
                evidence.snv_clinvar_evidence,
                evidence.sv_phenosv_evidence,
                evidence.sv_vcf,
                evidence.snv_vcf,
                evidence.hpo_path,
                evidence.age_of_onset
            )
        }
        ngs_prio(
            combined_priority_input,
            inheritance_mode,
            include_clinvar_report,
            allow_unphased_comphet
        )

        prio_report_input = ngs_prio.out.prio_vcf
            .join(ngs_prio.out.prio_gene_vcf, failOnMismatch: true, failOnDuplicate: true)
            .join(eh_filter.out, failOnMismatch: true, failOnDuplicate: true)
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
