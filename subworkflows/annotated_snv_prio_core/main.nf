include { multi_phen2gene } from '../../modules/multi_phen2gene/'
include { multi_rankscore } from '../../modules/multi_rankscore/'
include { multi_rankvar } from '../../modules/multi_rankvar/'
include { multi_phenotagger } from '../../modules/multi_phenotagger/'
include { multi_snp_prio } from '../../modules/multi_snp_prio/'
include { validate_preannotated_annovar_pair } from '../../modules/validate_preannotated_annovar_pair/'

// Shared core: pre-annotated ANNOVAR TXT/VCF SNP prioritization.
workflow ANNOTATED_SNV_PRIO_CORE {
	take:
	input_annotated_snv
	input_age
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

	main:
	validated_annovar = validate_preannotated_annovar_pair(
		input_annotated_snv.map { out_prefix, annovar_txt, annovar_vcf, phenotype_path, phenotype_format ->
			tuple(out_prefix, annovar_txt, annovar_vcf)
		}
	)

	clinical_note_input = input_annotated_snv
		.filter { out_prefix, annovar_txt, annovar_vcf, phenotype_path, phenotype_format -> phenotype_format == 'clinical_note' }
		.map { out_prefix, annovar_txt, annovar_vcf, phenotype_path, phenotype_format -> tuple(out_prefix, phenotype_path) }

	hpo_input = input_annotated_snv
		.filter { out_prefix, annovar_txt, annovar_vcf, phenotype_path, phenotype_format -> phenotype_format == 'hpo' }
		.map { out_prefix, annovar_txt, annovar_vcf, phenotype_path, phenotype_format -> tuple(out_prefix, phenotype_path) }

	phenotagger_result = multi_phenotagger(clinical_note_input)
	hpo_paths = phenotagger_result.mix(hpo_input)
	phen2gene_result = multi_phen2gene(hpo_paths)

	validated_annovar_txt = validated_annovar.map { out_prefix, annovar_txt, annovar_vcf -> tuple(out_prefix, annovar_txt) }
	join_annovar_phen2gene = validated_annovar_txt.join(phen2gene_result)
	join_annovar_hpo = join_annovar_phen2gene.join(hpo_paths)

	rankscore_result = multi_rankscore(join_annovar_phen2gene, gnomad, rankscore_filter, rankscore_softwares, gq, phen2gene_top_n)
	rankvar_result = multi_rankvar(join_annovar_hpo, gnomad, gq, ad, rankvar_filter)
	rankscore_rankvar_join = rankscore_result.join(rankvar_result)

	validated_annovar_vcf = validated_annovar.map { out_prefix, annovar_txt, annovar_vcf -> tuple(out_prefix, annovar_vcf) }
	snp_prio_input = rankscore_rankvar_join.join(validated_annovar_vcf)
	hpo_with_age = hpo_paths.join(input_age).map { out_prefix, hpo_path, age_of_onset -> tuple(out_prefix, hpo_path, age_of_onset) }
	snp_prio_input_hpo = snp_prio_input.join(hpo_with_age)
	multi_snp_prio(snp_prio_input_hpo, inheritance_mode, include_clinvar_report, allow_unphased_comphet)

	emit:
	validated_annovar = validated_annovar
	phen2gene = phen2gene_result
	rankscore = rankscore_result
	rankvar = rankvar_result
	prio_vcf = multi_snp_prio.out[0]
	prio_gene_vcf = multi_snp_prio.out[1]
}
