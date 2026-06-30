

include { sniffles } from '../../modules/sniffles/'
include { CNVpytor } from '../../modules/cnvpytor/'
include { merge_longread_sv_callers } from '../../modules/merge_longread_sv_callers/'
include { SURVIVOR } from '../../modules/survivor/'
include { PhenoSV } from '../../modules/phenosv/'
include { Phen2gene } from '../../modules/phen2gene/'
include { NanoRepeat } from '../../modules/nanorepeat/'
include { ANNOVAR_SV } from '../../modules/annovar_sv/'
include { common_sv_filter } from '../../modules/common_sv_filter/'
include { phenotagger } from '../../modules/phenotagger/'
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
		phenotagger(note,out_prefix)
		hpo=phenotagger.out
	}
	Phen2gene(hpo,out_prefix)
	sniffles(bam,out_prefix,ref_fa)
	cnvpytor_mode = params.cnvpytor ? params.cnvpytor.toString().trim().toLowerCase() : "no"
	if ( cnvpytor_mode == "yes" ) {
		cnvpytor_dummy_vcf = ref_fa.map { ref_tuple -> ref_tuple[0] }
		cnvpytor_reference_conf = Channel.value(params.cnvpytor_reference_conf ? file(params.cnvpytor_reference_conf) : [])
		CNVpytor(bam,out_prefix,ref_fa,Channel.value(params.cnvpytor_reference_genome),cnvpytor_reference_conf,cnvpytor_dummy_vcf,Channel.value("no"),Channel.value(params.cnvpytor_bin_sizes),Channel.value(params.cnvpytor_primary_bin),Channel.value(params.cnvpytor_min_size))
		sv_merge_inputs = sniffles.out.combine(CNVpytor.out.vcf).map { combined_vcfs ->
			combined_vcfs.flatten()
		}
		merge_longread_sv_callers(sv_merge_inputs,out_prefix)
		sv_vcf = merge_longread_sv_callers.out
	}
	else {
		sv_vcf = sniffles.out
	}
	ANNOVAR_SV(sv_vcf,out_prefix,Phen2gene.out,"null","called")
	annovar_sv_for_downstream = ANNOVAR_SV.out
	if ( params.common_sv_filter.toString().trim().toLowerCase() == "yes" ) {
		common_sv_filter(ANNOVAR_SV.out,out_prefix)
		annovar_sv_for_downstream = common_sv_filter.out.filtered_vcf
	}
	SURVIVOR(annovar_sv_for_downstream,out_prefix)
	PhenoSV(SURVIVOR.out,out_prefix,hpo)
	NanoRepeat(bam,out_prefix,ref_fa)
	sv_prio(out_prefix,PhenoSV.out,annovar_sv_for_downstream,hpo,inheritance_mode,include_clinvar_report,allow_unphased_comphet)
}
