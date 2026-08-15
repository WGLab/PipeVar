include { SURVIVOR } from '../../modules/survivor/'
include { PhenoSV } from '../../modules/phenosv/'
include { Phen2gene } from '../../modules/phen2gene/'
include { RankVar } from '../../modules/rankvar/'
include { ANNOVAR_SV } from '../../modules/annovar_sv/'
include { Rankscore_analysis_preannotated } from '../../modules/rankscore_analysis/'
include { phenotagger } from '../../modules/phenotagger/'
include { phenogpt2 } from '../../modules/phenogpt2/'
include { ngs_prio } from '../../modules/ngs_prio/'
include { common_sv_filter } from '../../modules/common_sv_filter/'
include { variant_html_report_no_repeat } from '../../modules/variant_html_report/'

// Single sample: imported ANNOVAR SNV and SV inputs without BAM/CRAM-dependent analysis.
workflow SINGLE_ANNOTATED_SNV_SV {
	take:
	annovar_txt
	snv_vcf
	annovar_sv_vcf
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
	hpo = phenotype
	if ( phenotype_format == 'clinical_note' || phenotype_format == 'yes' ) {
		if ( params.phenotype_extractor.toString().trim().toLowerCase() == "phenogpt2" ) {
			phenogpt2(phenotype, out_prefix)
			hpo = phenogpt2.out
		}
		else {
			phenotagger(phenotype, out_prefix)
			hpo = phenotagger.out
		}
	}
	Phen2gene(hpo, out_prefix)
	RankVar(annovar_txt, Phen2gene.out, hpo, out_prefix, gnomad, gq, ad, rankvar_filter)
	rankscore_result = Rankscore_analysis_preannotated(annovar_txt, snv_vcf, Phen2gene.out, out_prefix, gnomad, rankscore_filter, rankscore_softwares, gq, ad, phen2gene_top_n)

	ANNOVAR_SV(annovar_sv_vcf, out_prefix, "preannotated")
	annovar_sv_for_downstream = ANNOVAR_SV.out
	if ( params.common_sv_filter.toString().trim().toLowerCase() == "yes" ) {
		common_sv_filter(ANNOVAR_SV.out, out_prefix)
		annovar_sv_for_downstream = common_sv_filter.out.filtered_vcf
	}
	SURVIVOR(annovar_sv_for_downstream, out_prefix)
	PhenoSV(SURVIVOR.out.phenosv_inputs, out_prefix, hpo)
	ngs_prio(out_prefix, RankVar.out, rankscore_result.out.rankscore, rankscore_result.out.clinvar, PhenoSV.out, annovar_sv_for_downstream, snv_vcf, hpo, inheritance_mode, include_clinvar_report, allow_unphased_comphet)
	variant_html_report_no_repeat(out_prefix, ngs_prio.out[0], ngs_prio.out[1])

	emit:
	prio_vcf = ngs_prio.out[0]
	prio_gene_vcf = ngs_prio.out[1]
	html_report = variant_html_report_no_repeat.out
}
