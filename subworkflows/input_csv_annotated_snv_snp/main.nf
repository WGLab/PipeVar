include { ANNOTATED_SNV_PRIO_CORE } from '../annotated_snv_prio_core'

// CSV batch: pre-annotated ANNOVAR TXT/VCF SNP prioritization path.
workflow INPUT_CSV_ANNOTATED_SNV_SNP {
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
		denovo_filter,
		denovo_pedigree,
		denovo_role_column,
		denovo_family_column,
		denovo_vcf_sample_column,
		denovo_exclude_contigs,
		denovo_sv_min_reciprocal_overlap
	)

	emit:
	prio_vcf = ANNOTATED_SNV_PRIO_CORE.out.prio_vcf
	prio_gene_vcf = ANNOTATED_SNV_PRIO_CORE.out.prio_gene_vcf
}
