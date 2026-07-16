include { ANNOTATED_SNV_PRIO_CORE } from '../annotated_snv_prio_core'

// Single sample: pre-annotated ANNOVAR TXT/VCF SNP prioritization path.
workflow SINGLE_ANNOTATED_SNV_SNP {
	take:
	annovar_txt
	vcf
	out_prefix
	phenotype
	phenotype_format
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
	input_annotated_snv = out_prefix
		.combine(annovar_txt)
		.combine(vcf)
		.combine(phenotype)
		.combine(phenotype_format)
		.map { prefix, annovar_txt_file, annovar_vcf_file, phenotype_path, phenotype_kind ->
			tuple(prefix, annovar_txt_file, annovar_vcf_file, phenotype_path, phenotype_kind)
		}

	input_meta = out_prefix.map { prefix -> tuple(prefix, "", "unknown") }

	ANNOTATED_SNV_PRIO_CORE(
		input_annotated_snv,
		input_meta,
		rankscore_filter,
		rankscore_softwares,
		phen2gene_top_n,
		gnomad,
		gq,
		ad,
		rankvar_filter,
		inheritance_mode,
		include_clinvar_report,
		allow_unphased_comphet,
		"no",
		"null",
		"role",
		"family_id",
		"vcf_sample",
		"MT,M,chrM,chrMT",
		"0.50"
	)

	emit:
	prio_vcf = ANNOTATED_SNV_PRIO_CORE.out.prio_vcf
	prio_gene_vcf = ANNOTATED_SNV_PRIO_CORE.out.prio_gene_vcf
}
