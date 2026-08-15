include { multi_phen2gene } from '../../modules/multi_phen2gene/'
include { multi_rankscore_preannotated } from '../../modules/multi_rankscore/'
include { multi_rankvar } from '../../modules/multi_rankvar/'
include { multi_phenotagger } from '../../modules/multi_phenotagger/'
include { multi_phenogpt2 } from '../../modules/multi_phenogpt2/'
include { multi_ngs_prio } from '../../modules/multi_ngs_prio/'
include { multi_annovar_sv } from '../../modules/multi_annovar_sv/'
include { multi_common_sv_filter } from '../../modules/multi_common_sv_filter/'
include { multi_survivor } from '../../modules/multi_survivor/'
include { multi_phenosv } from '../../modules/multi_phenosv/'
include { multi_variant_html_report_no_repeat } from '../../modules/variant_html_report/'
include { DENOVO_SNV_FILTER_CORE } from '../denovo_snv_filter_core'
include { DENOVO_SV_FILTER_CORE } from '../denovo_sv_filter_core'

// CSV batch: imported ANNOVAR SNV and SV inputs without BAM/CRAM-dependent analysis.
workflow INPUT_CSV_ANNOTATED_SNV_SV {
	take:
	input_annotated_snv_sv
	input_meta
	rankscore_filter
	rankscore_softwares
	phen2gene_top_n
	gnomad
	gq
	ad
	rankvar_filter
	inheritance_mode
	include_clinvar_report
	allow_unphased_comphet
	denovo_filter
	denovo_pedigree
	denovo_role_column
	denovo_family_column
	denovo_vcf_sample_column
	denovo_exclude_contigs
	denovo_sv_min_reciprocal_overlap

	main:
	annovar_for_downstream = input_annotated_snv_sv.map { out_prefix, annovar_txt, annovar_vcf, annovar_sv_vcf, phenotype_path, phenotype_format ->
		tuple(out_prefix, annovar_txt, annovar_vcf)
	}
	analysis_input = input_annotated_snv_sv
	if ( denovo_filter == "yes" ) {
		denovo_snv_result = DENOVO_SNV_FILTER_CORE(annovar_for_downstream, denovo_pedigree, denovo_role_column, denovo_family_column, denovo_vcf_sample_column, denovo_exclude_contigs)
		annovar_for_downstream = denovo_snv_result.records
		proband_keys = denovo_snv_result.records.map { out_prefix, txt_file, vcf_file -> tuple(out_prefix, true) }
		analysis_input = input_annotated_snv_sv.join(proband_keys, failOnDuplicate: true).map { out_prefix, annovar_txt, annovar_vcf, annovar_sv_vcf, phenotype_path, phenotype_format, marker ->
			tuple(out_prefix, annovar_txt, annovar_vcf, annovar_sv_vcf, phenotype_path, phenotype_format)
		}
	}

	clinical_note_input = analysis_input
		.filter { out_prefix, annovar_txt, annovar_vcf, annovar_sv_vcf, phenotype_path, phenotype_format -> phenotype_format == 'clinical_note' }
		.map { out_prefix, annovar_txt, annovar_vcf, annovar_sv_vcf, phenotype_path, phenotype_format -> tuple(out_prefix, phenotype_path) }
	hpo_input = analysis_input
		.filter { out_prefix, annovar_txt, annovar_vcf, annovar_sv_vcf, phenotype_path, phenotype_format -> phenotype_format == 'hpo' }
		.map { out_prefix, annovar_txt, annovar_vcf, annovar_sv_vcf, phenotype_path, phenotype_format -> tuple(out_prefix, phenotype_path) }
	if ( params.phenotype_extractor.toString().trim().toLowerCase() == "phenogpt2" ) {
		phenotype_extractor_result = multi_phenogpt2(clinical_note_input)
	}
	else {
		phenotype_extractor_result = multi_phenotagger(clinical_note_input)
	}
	hpo_paths = phenotype_extractor_result.mix(hpo_input)
	phen2gene_result = multi_phen2gene(hpo_paths)

	annovar_txt_for_rank = annovar_for_downstream.map { out_prefix, annovar_txt, annovar_vcf -> tuple(out_prefix, annovar_txt) }
	join_annovar_phen2gene = annovar_txt_for_rank.join(phen2gene_result, failOnMismatch: true, failOnDuplicate: true)
	join_annovar_hpo = join_annovar_phen2gene.join(hpo_paths, failOnMismatch: true, failOnDuplicate: true)
	annovar_vcf_for_rankscore = annovar_for_downstream.map { out_prefix, annovar_txt, annovar_vcf -> tuple(out_prefix, annovar_vcf) }
	rankscore_input = join_annovar_phen2gene
		.join(annovar_vcf_for_rankscore, failOnMismatch: true, failOnDuplicate: true)
		.map { out_prefix, annovar_txt, phen2gene, annovar_vcf -> tuple(out_prefix, annovar_txt, annovar_vcf, phen2gene) }
	rankscore_result = multi_rankscore_preannotated(rankscore_input, gnomad, rankscore_filter, rankscore_softwares, gq, ad, phen2gene_top_n)
	rankvar_result = multi_rankvar(join_annovar_hpo, gnomad, gq, ad, rankvar_filter)

	annotated_sv_input = input_annotated_snv_sv.map { out_prefix, annovar_txt, annovar_vcf, annovar_sv_vcf, phenotype_path, phenotype_format ->
		tuple(out_prefix, annovar_sv_vcf)
	}
	sv_for_annotation = annotated_sv_input
	if ( denovo_filter == "yes" ) {
		denovo_sv_result = DENOVO_SV_FILTER_CORE(annotated_sv_input, denovo_pedigree, denovo_role_column, denovo_family_column, denovo_vcf_sample_column, denovo_exclude_contigs, denovo_sv_min_reciprocal_overlap)
		sv_for_annotation = denovo_sv_result.records
	}
	sv_result_annovar = sv_for_annotation.map { out_prefix, vcf_file ->
		tuple(out_prefix, vcf_file, "preannotated")
	}
	annovar_sv_result = multi_annovar_sv(sv_result_annovar)
	annovar_sv_for_downstream = annovar_sv_result
	if ( params.common_sv_filter.toString().trim().toLowerCase() == "yes" ) {
		multi_common_sv_filter(annovar_sv_for_downstream)
		annovar_sv_for_downstream = multi_common_sv_filter.out.filtered_vcf
	}
	survivor_result = multi_survivor(annovar_sv_for_downstream)
	phenosv_input = survivor_result.join(hpo_paths, failOnMismatch: true, failOnDuplicate: true)
	phenosv_result = multi_phenosv(phenosv_input)

	annovar_vcf_for_prio = annovar_for_downstream.map { out_prefix, annovar_txt, annovar_vcf -> tuple(out_prefix, annovar_vcf) }
	phenosv_annovar_snv = phenosv_result.join(annovar_vcf_for_prio, failOnMismatch: true, failOnDuplicate: true)
	sv_join = phenosv_annovar_snv.join(annovar_sv_for_downstream, failOnMismatch: true, failOnDuplicate: true)
	rankscore_join = sv_join.join(rankscore_result, failOnMismatch: true, failOnDuplicate: true)
	rankvar_join = rankscore_join.join(rankvar_result, failOnMismatch: true, failOnDuplicate: true)
	hpo_with_age = hpo_paths.join(input_meta, failOnMismatch: true, failOnDuplicate: true).map { out_prefix, hpo_path, age_of_onset, sex -> tuple(out_prefix, hpo_path, age_of_onset, sex) }
	rankvar_join_hpo = rankvar_join.join(hpo_with_age, failOnMismatch: true, failOnDuplicate: true)
	rankvar_join_hpo_ordered = rankvar_join_hpo.map { out_prefix, sv_pathogenic, snv_vcf_path, sv_vcf_path, snv_rankscore, snv_pathogenic, snv_rankvar, hpo_path, age_of_onset, sex ->
		tuple(out_prefix, snv_rankvar, snv_rankscore, snv_pathogenic, sv_pathogenic, sv_vcf_path, snv_vcf_path, hpo_path, age_of_onset, sex)
	}
	multi_ngs_prio(rankvar_join_hpo_ordered, inheritance_mode, include_clinvar_report, allow_unphased_comphet)

	prio_report_input = multi_ngs_prio.out[0]
		.map { prio_vcf -> tuple(prio_vcf.name.replaceFirst(/\.prio\.vcf$/, ""), prio_vcf) }
		.join(
			multi_ngs_prio.out[1].map { prio_gene_report ->
				tuple(prio_gene_report.name.replaceFirst(/\.prio_gene\.vcf$/, ""), prio_gene_report)
			}
		)
	multi_variant_html_report_no_repeat(prio_report_input)

	emit:
	prio_vcf = multi_ngs_prio.out[0]
	prio_gene_vcf = multi_ngs_prio.out[1]
	html_report = multi_variant_html_report_no_repeat.out
}
