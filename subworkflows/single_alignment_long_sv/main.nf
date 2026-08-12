

include { sniffles } from '../../modules/sniffles/'
include { SURVIVOR } from '../../modules/survivor/'
include { PhenoSV } from '../../modules/phenosv/'
include { NanoRepeat } from '../../modules/nanorepeat/'
include { ANNOVAR_SV } from '../../modules/annovar_sv/'
include { common_sv_filter } from '../../modules/common_sv_filter/'
include { phenotagger } from '../../modules/phenotagger/'
include { phenogpt2 } from '../../modules/phenogpt2/'
include { sv_prio } from '../../modules/sv_prio/'


// Single sample: long-read SV/STR prioritization path with Sniffles -> ANNOVAR_SV.
workflow SINGLE_ALIGNMENT_LONG_SV {
	take:
	bam
	out_prefix
	ref_fa
	note
	is_note
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
	sniffles(bam,out_prefix,ref_fa)
	ANNOVAR_SV(sniffles.out,out_prefix,"called")
	annovar_sv_for_downstream = ANNOVAR_SV.out
	if ( params.common_sv_filter.toString().trim().toLowerCase() == "yes" ) {
		common_sv_filter(ANNOVAR_SV.out,out_prefix)
		annovar_sv_for_downstream = common_sv_filter.out.filtered_vcf
	}
	SURVIVOR(annovar_sv_for_downstream,out_prefix)
	PhenoSV(SURVIVOR.out.phenosv_inputs,out_prefix,hpo)
	NanoRepeat(bam,out_prefix,ref_fa)
	sv_prio(out_prefix,PhenoSV.out,annovar_sv_for_downstream,hpo,inheritance_mode,include_clinvar_report,allow_unphased_comphet)
}
