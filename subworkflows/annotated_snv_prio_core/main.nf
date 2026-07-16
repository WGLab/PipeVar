include { multi_phen2gene } from '../../modules/multi_phen2gene/'
include { multi_rankscore } from '../../modules/multi_rankscore/'
include { multi_rankvar } from '../../modules/multi_rankvar/'
include { multi_phenotagger } from '../../modules/multi_phenotagger/'
include { multi_phenogpt2 } from '../../modules/multi_phenogpt2/'
include { multi_snp_prio } from '../../modules/multi_snp_prio/'
include { validate_preannotated_annovar_pair } from '../../modules/validate_preannotated_annovar_pair/'
include { DENOVO_SNV_FILTER_CORE } from '../denovo_snv_filter_core'

// Shared core: pre-annotated ANNOVAR TXT/VCF SNP prioritization.
workflow ANNOTATED_SNV_PRIO_CORE {
	take:
	input_annotated_snv
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
	validated_annovar = validate_preannotated_annovar_pair(
		input_annotated_snv.map { out_prefix, annovar_txt, annovar_vcf, phenotype_path, phenotype_format ->
			tuple(out_prefix, annovar_txt, annovar_vcf)
		}
	)
	validated_for_downstream = validated_annovar
	analysis_input = input_annotated_snv
	if ( denovo_filter == "yes" ) {
		denovo_result = DENOVO_SNV_FILTER_CORE(validated_annovar, denovo_pedigree, denovo_role_column, denovo_family_column, denovo_vcf_sample_column, denovo_exclude_contigs)
		validated_for_downstream = denovo_result.records
		proband_keys = denovo_result.records.map { out_prefix, txt_file, vcf_file -> tuple(out_prefix, true) }
		analysis_input = input_annotated_snv.join(proband_keys, failOnDuplicate: true).map { out_prefix, annovar_txt, annovar_vcf, phenotype_path, phenotype_format, marker ->
			tuple(out_prefix, annovar_txt, annovar_vcf, phenotype_path, phenotype_format)
		}
	}

	clinical_note_input = analysis_input
		.filter { out_prefix, annovar_txt, annovar_vcf, phenotype_path, phenotype_format -> phenotype_format == 'clinical_note' }
		.map { out_prefix, annovar_txt, annovar_vcf, phenotype_path, phenotype_format -> tuple(out_prefix, phenotype_path) }

	hpo_input = analysis_input
		.filter { out_prefix, annovar_txt, annovar_vcf, phenotype_path, phenotype_format -> phenotype_format == 'hpo' }
		.map { out_prefix, annovar_txt, annovar_vcf, phenotype_path, phenotype_format -> tuple(out_prefix, phenotype_path) }

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

	rankscore_result = multi_rankscore(join_annovar_phen2gene, gnomad, rankscore_filter, rankscore_softwares, gq, phen2gene_top_n)
	rankvar_result = multi_rankvar(join_annovar_hpo, gnomad, gq, ad, rankvar_filter)
	rankscore_rankvar_join = rankscore_result.join(rankvar_result)

	validated_annovar_vcf = validated_for_downstream.map { out_prefix, annovar_txt, annovar_vcf -> tuple(out_prefix, annovar_vcf) }
	snp_prio_input = rankscore_rankvar_join.join(validated_annovar_vcf)
	hpo_with_age = hpo_paths.join(input_meta, failOnMismatch: true, failOnDuplicate: true).map { out_prefix, hpo_path, age_of_onset, sex -> tuple(out_prefix, hpo_path, age_of_onset, sex) }
	snp_prio_input_hpo = snp_prio_input.join(hpo_with_age)
	multi_snp_prio(snp_prio_input_hpo, inheritance_mode, include_clinvar_report, allow_unphased_comphet)

	emit:
	validated_annovar = validated_for_downstream
	phen2gene = phen2gene_result
	rankscore = rankscore_result
	rankvar = rankvar_result
	prio_vcf = multi_snp_prio.out[0]
	prio_gene_vcf = multi_snp_prio.out[1]
}
