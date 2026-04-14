

include { SURVIVOR } from '../../modules/survivor/'
include { PhenoSV } from '../../modules/phenosv/'
include { Phen2gene } from '../../modules/phen2gene/'
include { Manta } from '../../modules/manta/'
include { MELT } from '../../modules/melt/'
include { normalize_shortread_alignment } from '../../modules/normalize_shortread_alignment/'
include { CNVnator } from '../../modules/cnvnator/'
include { merge_shortread_sv_callers } from '../../modules/merge_shortread_sv_callers/'
include { ANNOVAR_SV } from '../../modules/annovar_sv/'
include { ExpansionHunter } from '../../modules/expansion_hunter/'
include { phenotagger } from '../../modules/phenotagger/'
include { eh_filter } from '../../modules/eh_filter/'
include { sv_prio } from '../../modules/sv_prio/'


// Single sample: short-read SV/STR prioritization path with Manta -> ANNOVAR_SV.
workflow SINGLE_ALIGNMENT_NGS_SV {
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
	Manta(bam,out_prefix,ref_fa)
	melt_mode = params.melt ? params.melt.toString().trim().toLowerCase() : "no"
	if ( melt_mode == "yes" ) {
		MELT(bam,out_prefix,ref_fa)
	}

	cnvnator_mode = params.cnvnator ? params.cnvnator.toString().trim().toLowerCase() : "yes"
	if ( cnvnator_mode != "no" && melt_mode == "yes" ) {
		normalize_shortread_alignment(bam,out_prefix,ref_fa)
		CNVnator(normalize_shortread_alignment.out,out_prefix,ref_fa,params.cnvnator_bin_size)
		sv_merge_inputs = Manta.out.combine(CNVnator.out.vcf).combine(MELT.out).map { combined_vcfs ->
			combined_vcfs.flatten()
		}
		merge_shortread_sv_callers(sv_merge_inputs,out_prefix)
		sv_vcf=merge_shortread_sv_callers.out
	}
	else if ( cnvnator_mode != "no" ) {
		normalize_shortread_alignment(bam,out_prefix,ref_fa)
		CNVnator(normalize_shortread_alignment.out,out_prefix,ref_fa,params.cnvnator_bin_size)
		sv_merge_inputs = Manta.out.combine(CNVnator.out.vcf).map { combined_vcfs ->
			combined_vcfs.flatten()
		}
		merge_shortread_sv_callers(sv_merge_inputs,out_prefix)
		sv_vcf=merge_shortread_sv_callers.out
	}
	else if ( melt_mode == "yes" ) {
		sv_merge_inputs = Manta.out.combine(MELT.out).map { combined_vcfs ->
			combined_vcfs.flatten()
		}
		merge_shortread_sv_callers(sv_merge_inputs,out_prefix)
		sv_vcf=merge_shortread_sv_callers.out
	}
	else {
		sv_vcf=Manta.out
	}

	ANNOVAR_SV(sv_vcf,out_prefix,Phen2gene.out,"null")
	SURVIVOR(ANNOVAR_SV.out,out_prefix)
	PhenoSV(SURVIVOR.out,out_prefix,hpo)
	ExpansionHunter(bam,out_prefix,ref_fa)
	eh_filter(out_prefix,ExpansionHunter.out)
		sv_prio(out_prefix,PhenoSV.out,ANNOVAR_SV.out,hpo,inheritance_mode,include_clinvar_report,allow_unphased_comphet)
}
