

include { SURVIVOR } from '../../modules/survivor/'
include { PhenoSV } from '../../modules/phenosv/'
include { Phen2gene } from '../../modules/phen2gene/'
include { ANNOVAR_SV } from '../../modules/annovar_sv/'
include { common_sv_filter } from '../../modules/common_sv_filter/'
include { phenotagger } from '../../modules/phenotagger/'
include { phenogpt2 } from '../../modules/phenogpt2/'
include { phen2gene_filter } from '../../modules/reduce_region_phen2gene/'
include { sv_prio } from '../../modules/sv_prio/'


// Single sample: VCF SV re-annotation and SV prioritization path.
workflow SINGLE_ALIGNMENT_VCF_SV {
	take:
	vcf
	out_prefix
	ref_fa
	note
	phen2gene_top_n
	is_note
	target
	inheritance_mode
	include_clinvar_report
	allow_unphased_comphet

	main:
	
		hpo=note
		if ( is_note == "yes" ) {
		if ( params.phenotype_extractor.toString().trim().toLowerCase() == "phenogpt2" ) {
			phenogpt2(note,out_prefix)
			hpo=phenogpt2.out
		}
		else {
			phenotagger(note,out_prefix)
			hpo=phenotagger.out
		}
	}
		Phen2gene(hpo,out_prefix)
		if ( target == "yes" ) {
			phen2_gene_bed=phen2gene_filter(Phen2gene.out,ref_fa,out_prefix,phen2gene_top_n)
			annovar_sv_bed = phen2_gene_bed
		}
		else {
			annovar_sv_bed = target
		}
		ANNOVAR_SV(vcf,out_prefix,Phen2gene.out,annovar_sv_bed,"called")
		annovar_sv_for_downstream = ANNOVAR_SV.out
		if ( params.common_sv_filter.toString().trim().toLowerCase() == "yes" ) {
			common_sv_filter(ANNOVAR_SV.out,out_prefix)
			annovar_sv_for_downstream = common_sv_filter.out.filtered_vcf
		}
		SURVIVOR(annovar_sv_for_downstream,out_prefix)
		PhenoSV(SURVIVOR.out.phenosv_inputs,out_prefix,hpo)
		sv_prio(out_prefix,PhenoSV.out,annovar_sv_for_downstream,hpo,inheritance_mode,include_clinvar_report,allow_unphased_comphet)
}
