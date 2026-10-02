include { annovar_sv } from '../../modules/annovar_sv/'
include { common_sv_filter } from '../../modules/common_sv_filter/'
include { phen2gene } from '../../modules/phen2gene/'
include { phenogpt2 } from '../../modules/phenogpt2/'
include { phenosv } from '../../modules/phenosv/'
include { phenotagger } from '../../modules/phenotagger/'
include { rankscore_preannotated } from '../../modules/rankscore_analysis/'
include { rankvar } from '../../modules/rankvar/'
include { survivor } from '../../modules/survivor/'

// Shared validation, phenotype normalization, and ranking for imported ANNOVAR
// SNV TXT/VCF pairs. Route wrappers adapt their distinct manifests to this
// canonical five-field record before invoking the helper.
workflow PREANNOTATED_SNV_EVIDENCE {
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

    main:
    annovar_for_downstream = input_annotated_snv.map { sample_id, annovar_txt, annovar_vcf, phenotype_path, phenotype_format ->
        tuple(sample_id, annovar_txt, annovar_vcf)
    }
    clinical_note_input = input_annotated_snv
        .filter { sample_id, annovar_txt, annovar_vcf, phenotype_path, phenotype_format -> phenotype_format == 'clinical_note' }
        .map { sample_id, annovar_txt, annovar_vcf, phenotype_path, phenotype_format -> tuple(sample_id, phenotype_path) }
    hpo_input = input_annotated_snv
        .filter { sample_id, annovar_txt, annovar_vcf, phenotype_path, phenotype_format -> phenotype_format == 'hpo' }
        .map { sample_id, annovar_txt, annovar_vcf, phenotype_path, phenotype_format -> tuple(sample_id, phenotype_path) }

    if (params.phenotype_extractor.toString().trim().toLowerCase() == 'phenogpt2') {
        phenotype_extractor_result = phenogpt2(clinical_note_input)
    }
    else {
        phenotype_extractor_result = phenotagger(clinical_note_input)
    }
    hpo_paths = phenotype_extractor_result.mix(hpo_input)
    phen2gene_result = phen2gene(hpo_paths)

    annovar_txt_for_rank = annovar_for_downstream.map { sample_id, annovar_txt, annovar_vcf -> tuple(sample_id, annovar_txt) }
    annovar_with_ranked_genes = annovar_txt_for_rank.join(
        phen2gene_result,
        failOnMismatch: true,
        failOnDuplicate: true
    )
    annovar_with_ranked_genes_and_hpo = annovar_with_ranked_genes.join(
        hpo_paths,
        failOnMismatch: true,
        failOnDuplicate: true
    )
    annovar_vcf_for_rankscore = annovar_for_downstream.map { sample_id, annovar_txt, annovar_vcf -> tuple(sample_id, annovar_vcf) }
    rankscore_input = annovar_with_ranked_genes
        .join(annovar_vcf_for_rankscore, failOnMismatch: true, failOnDuplicate: true)
        .map { sample_id, annovar_txt, ranked_genes, annovar_vcf -> tuple(sample_id, annovar_txt, annovar_vcf, ranked_genes) }
    rankscore_result = rankscore_preannotated(
        rankscore_input,
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
    phenotype_metadata = hpo_paths
        .join(clinical_metadata_by_sample, failOnMismatch: true, failOnDuplicate: true)
        .map { sample_id, hpo_path, age_of_onset -> tuple(sample_id, hpo_path, age_of_onset) }

    emit:
    annovar = annovar_for_downstream
    hpo = hpo_paths
    phenotype_metadata = phenotype_metadata
    phen2gene = phen2gene_result
    rankscore = rankscore_result
    rankvar = rankvar_result
}

// Shared imported-SV annotation and phenotype scoring. Called-SV routes keep
// their separate caller/merge chain and use annovar_sv in called mode.
workflow PREANNOTATED_IMPORTED_SV_CHAIN {
    take:
    imported_sv_by_sample
    hpo_files_by_sample

    main:
    sv_result_annovar = imported_sv_by_sample.map { sample_id, vcf_file ->
        tuple(sample_id, vcf_file, 'preannotated')
    }
    annovar_sv_for_downstream = annovar_sv(sv_result_annovar)
    if (params.common_sv_filter.toString().trim().toLowerCase() == 'yes') {
        common_sv_filter(annovar_sv_for_downstream)
        annovar_sv_for_downstream = common_sv_filter.out.filtered_vcf
    }
    survivor_result = survivor(annovar_sv_for_downstream)
    phenosv_input = survivor_result.join(hpo_files_by_sample, failOnMismatch: true, failOnDuplicate: true)
    phenosv_result = phenosv(phenosv_input)

    emit:
    annovar_sv = annovar_sv_for_downstream
    phenosv = phenosv_result
}

// Assemble the common nine-field ngs_prio record from independently produced
// SNV, SV, phenotype, and clinical evidence. Each join remains strict and
// sample-keyed so missing or duplicate evidence fails at its producer boundary.
workflow PREANNOTATED_COMBINED_EVIDENCE {
    take:
    annovar_snv
    rankscore_result
    rankvar_result
    phenosv_result
    annovar_sv
    phenotype_metadata

    main:
    annovar_vcf_for_prio = annovar_snv.map { sample_id, annovar_txt, annovar_vcf ->
        tuple(sample_id, annovar_vcf)
    }
    phenosv_annovar_snv = phenosv_result.join(
        annovar_vcf_for_prio,
        failOnMismatch: true,
        failOnDuplicate: true
    )
    sv_join = phenosv_annovar_snv.join(annovar_sv, failOnMismatch: true, failOnDuplicate: true)
    rankscore_join = sv_join.join(rankscore_result, failOnMismatch: true, failOnDuplicate: true)
    rankvar_join = rankscore_join.join(rankvar_result, failOnMismatch: true, failOnDuplicate: true)
    combined_evidence_by_sample = rankvar_join
        .join(phenotype_metadata, failOnMismatch: true, failOnDuplicate: true)
        .map { sample_id, sv_phenosv_evidence, snv_vcf, sv_vcf, snv_rankscore, snv_clinvar_evidence, snv_rankvar, hpo_path, age_of_onset ->
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
    priority_input = combined_evidence_by_sample.map { sample_id, evidence ->
        tuple(sample_id, evidence.snv_rankvar, evidence.snv_rankscore, evidence.snv_clinvar_evidence,
            evidence.sv_phenosv_evidence, evidence.sv_vcf, evidence.snv_vcf, evidence.hpo_path,
            evidence.age_of_onset)
    }

    emit:
    priority_input
}
