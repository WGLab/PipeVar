include { multi_phen2gene } from '../../modules/multi_phen2gene/'
include { multi_rankscore } from '../../modules/multi_rankscore/'
include { multi_rankvar } from '../../modules/multi_rankvar/'
include { multi_expansionhunter } from '../../modules/multi_expansionhunter/'
include { multi_eh_filter } from '../../modules/multi_eh_filter/'
include { multi_phenotagger } from '../../modules/multi_phenotagger/'
include { multi_phenogpt2 } from '../../modules/multi_phenogpt2/'
include { multi_ngs_prio } from '../../modules/multi_ngs_prio/'
include { multi_annovar_sv } from '../../modules/multi_annovar_sv/'
include { multi_common_sv_filter } from '../../modules/multi_common_sv_filter/'
include { multi_survivor } from '../../modules/multi_survivor/'
include { multi_phenosv } from '../../modules/multi_phenosv/'
include { validate_preannotated_annovar_pair } from '../../modules/validate_preannotated_annovar_pair/'
include { multi_mito_prep_mutect2 } from '../../modules/multi_mito_prep_mutect2/'
include { multi_mito_mutect2 } from '../../modules/multi_mito_mutect2/'
include { multi_mito_annotation } from '../../modules/multi_mito_annotation/'
include { multi_mito_prio } from '../../modules/multi_mito_prio/'
include { multi_variant_html_report; multi_variant_html_report_with_mito } from '../../modules/variant_html_report/'
include { DENOVO_SNV_FILTER_CORE } from '../denovo_snv_filter_core'
include { DENOVO_SV_FILTER_CORE } from '../denovo_sv_filter_core'

// CSV batch: imported ANNOVAR SNV + short-read SV/CNV/(optional) mito analysis.
workflow INPUT_CSV_ANNOTATED_ALL_NGS {
	take:
	input_annotated_ngs
	input_meta
	ref_fa
	eh_ref_fa
	eh_variant_catalog
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
	mito_ref_fa
	mito_contig
	denovo_filter
	denovo_pedigree
	denovo_role_column
	denovo_family_column
	denovo_vcf_sample_column
	denovo_exclude_contigs
	denovo_sv_min_reciprocal_overlap

	main:
	validate_input = input_annotated_ngs.map { out_prefix, annovar_txt, annovar_vcf, annovar_sv_vcf, bam_file, bai_file, phenotype_path, phenotype_format ->
		tuple(out_prefix, annovar_txt, annovar_vcf)
	}
	validated_annovar = validate_preannotated_annovar_pair(validate_input)
	validated_for_downstream = validated_annovar
	analysis_input = input_annotated_ngs
	if ( denovo_filter == "yes" ) {
		denovo_snv_result = DENOVO_SNV_FILTER_CORE(validated_annovar, denovo_pedigree, denovo_role_column, denovo_family_column, denovo_vcf_sample_column, denovo_exclude_contigs)
		validated_for_downstream = denovo_snv_result.records
		proband_keys = denovo_snv_result.records.map { out_prefix, txt_file, vcf_file -> tuple(out_prefix, true) }
		analysis_input = input_annotated_ngs.join(proband_keys, failOnDuplicate: true).map { out_prefix, annovar_txt, annovar_vcf, annovar_sv_vcf, bam_file, bai_file, phenotype_path, phenotype_format, marker ->
			tuple(out_prefix, annovar_txt, annovar_vcf, annovar_sv_vcf, bam_file, bai_file, phenotype_path, phenotype_format)
		}
	}

	clinical_note_input = analysis_input
		.filter { out_prefix, annovar_txt, annovar_vcf, annovar_sv_vcf, bam_file, bai_file, phenotype_path, phenotype_format -> phenotype_format == 'clinical_note' }
		.map { out_prefix, annovar_txt, annovar_vcf, annovar_sv_vcf, bam_file, bai_file, phenotype_path, phenotype_format -> tuple(out_prefix, phenotype_path) }
	hpo_input = analysis_input
		.filter { out_prefix, annovar_txt, annovar_vcf, annovar_sv_vcf, bam_file, bai_file, phenotype_path, phenotype_format -> phenotype_format == 'hpo' }
		.map { out_prefix, annovar_txt, annovar_vcf, annovar_sv_vcf, bam_file, bai_file, phenotype_path, phenotype_format -> tuple(out_prefix, phenotype_path) }
	if ( params.phenotype_extractor.toString().trim().toLowerCase() == "phenogpt2" ) {
		phenotype_extractor_result = multi_phenogpt2(clinical_note_input)
	}
	else {
		phenotype_extractor_result = multi_phenotagger(clinical_note_input)
	}
	hpo_paths = phenotype_extractor_result.mix(hpo_input)
	phen2gene_result = multi_phen2gene(hpo_paths)

	validated_annovar_txt = validated_for_downstream.map { out_prefix, annovar_txt, annovar_vcf -> tuple(out_prefix, annovar_txt) }
	join_annovar_phen2gene = validated_annovar_txt.join(phen2gene_result)
	join_annovar_hpo = join_annovar_phen2gene.join(hpo_paths)
	rankscore_result = multi_rankscore(join_annovar_phen2gene, gnomad, rankscore_filter, rankscore_softwares, gq, ad, phen2gene_top_n)
	rankvar_result = multi_rankvar(join_annovar_hpo, gnomad, gq, ad, rankvar_filter)

	input_bam_with_bam = analysis_input.map { out_prefix, annovar_txt, annovar_vcf, annovar_sv_vcf, bam_file, bai_file, phenotype_path, phenotype_format ->
		tuple(out_prefix, bam_file, bai_file)
	}
	multi_eh_result = multi_expansionhunter(input_bam_with_bam, eh_ref_fa, eh_variant_catalog)
	multi_eh_filter(multi_eh_result)

	annotated_sv_input = input_annotated_ngs.map { out_prefix, annovar_txt, annovar_vcf, annovar_sv_vcf, bam_file, bai_file, phenotype_path, phenotype_format ->
		tuple(out_prefix, annovar_sv_vcf)
	}
	sv_for_annotation = annotated_sv_input
	if ( denovo_filter == "yes" ) {
		denovo_sv_result = DENOVO_SV_FILTER_CORE(annotated_sv_input, denovo_pedigree, denovo_role_column, denovo_family_column, denovo_vcf_sample_column, denovo_exclude_contigs, denovo_sv_min_reciprocal_overlap)
		sv_for_annotation = denovo_sv_result.records
	}
	sv_result_annovar = sv_for_annotation.join(phen2gene_result, failOnMismatch: true, failOnDuplicate: true).map { out_prefix, vcf_file, phen2gene_file ->
		tuple(out_prefix, vcf_file, phen2gene_file, "null", "preannotated")
	}
	annovar_sv_result = multi_annovar_sv(sv_result_annovar)
	annovar_sv_for_downstream = annovar_sv_result
	if ( params.common_sv_filter.toString().trim().toLowerCase() == "yes" ) {
		multi_common_sv_filter(annovar_sv_for_downstream)
		annovar_sv_for_downstream = multi_common_sv_filter.out.filtered_vcf
	}
	survivor_result = multi_survivor(annovar_sv_for_downstream)
	phenosv_input = survivor_result.join(hpo_paths)
	phenosv_result = multi_phenosv(phenosv_input)

	validated_annovar_vcf = validated_for_downstream.map { out_prefix, annovar_txt, annovar_vcf -> tuple(out_prefix, annovar_vcf) }
	phenosv_annovar_snv = phenosv_result.join(validated_annovar_vcf)
	sv_join = phenosv_annovar_snv.join(annovar_sv_for_downstream)
	rankscore_join = sv_join.join(rankscore_result)
	rankvar_join = rankscore_join.join(rankvar_result)
	hpo_with_age = hpo_paths.join(input_meta, failOnMismatch: true, failOnDuplicate: true).map { out_prefix, hpo_path, age_of_onset, sex -> tuple(out_prefix, hpo_path, age_of_onset, sex) }
	rankvar_join_hpo = rankvar_join.join(hpo_with_age)
	rankvar_join_hpo_ordered = rankvar_join_hpo.map { out_prefix, sv_pathogenic, snv_vcf_path, sv_vcf_path, snv_rankscore, snv_pathogenic, snv_rankvar, hpo_path, age_of_onset, sex ->
		tuple(out_prefix, snv_rankvar, snv_rankscore, snv_pathogenic, sv_pathogenic, sv_vcf_path, snv_vcf_path, hpo_path, age_of_onset, sex)
	}
	multi_ngs_prio(rankvar_join_hpo_ordered, inheritance_mode, include_clinvar_report, allow_unphased_comphet)

	if ( mito_ref_fa != null ) {
		mutect2_ref_fa = mito_ref_fa.map { fa_file, fai_file, dict_file, bwa_amb, bwa_ann, bwa_bwt, bwa_pac, bwa_sa ->
			tuple(fa_file, fai_file, dict_file)
		}
		prepped = multi_mito_prep_mutect2(input_bam_with_bam, mito_ref_fa, mito_contig)
		mt_vcf = multi_mito_mutect2(prepped, mutect2_ref_fa, mito_contig)
		annotated = multi_mito_annotation(mt_vcf)
		multi_mito_prio(annotated)
	}

	prio_report_input = multi_ngs_prio.out[0]
		.map { prio_vcf ->
			def prefix = prio_vcf.name.replaceFirst(/\.prio\.vcf$/, "")
			tuple(prefix, prio_vcf)
		}
		.join(
			multi_ngs_prio.out[1].map { prio_gene_report ->
				def prefix = prio_gene_report.name.replaceFirst(/\.prio_gene\.vcf$/, "")
				tuple(prefix, prio_gene_report)
			}
		)
		.join(
			multi_eh_filter.out.map { repeat_tsv ->
				def prefix = repeat_tsv.name.replaceFirst(/\.eh\.tsv$/, "")
				tuple(prefix, repeat_tsv)
			}
		)
		.map { out_prefix, prio_vcf, prio_gene_report, repeat_tsv ->
			tuple(out_prefix, prio_vcf, prio_gene_report, repeat_tsv)
		}

	if ( mito_ref_fa != null ) {
		prio_report_input_with_mito = prio_report_input
			.join(multi_mito_prio.out)
			.map { out_prefix, prio_vcf, prio_gene_report, repeat_tsv, mito_tsv ->
				tuple(out_prefix, prio_vcf, prio_gene_report, repeat_tsv, mito_tsv)
			}
		multi_variant_html_report_with_mito(prio_report_input_with_mito)
	}
	else {
		multi_variant_html_report(prio_report_input)
	}

	emit:
	prio_vcf = multi_ngs_prio.out[0]
	prio_gene_vcf = multi_ngs_prio.out[1]
	mito_vcf = mito_ref_fa != null ? mt_vcf : Channel.empty()
	mito_annotated_tsv = mito_ref_fa != null ? annotated : Channel.empty()
	mito_prioritized_tsv = mito_ref_fa != null ? multi_mito_prio.out : Channel.empty()
}
