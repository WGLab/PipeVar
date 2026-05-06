

include { SURVIVOR } from '../../modules/survivor/'
include { PhenoSV } from '../../modules/phenosv/'
include { Phen2gene } from '../../modules/phen2gene/'
include { Manta } from '../../modules/manta/'
include { scramble } from '../../modules/scramble/'
include { scramble_ref_prep } from '../../modules/scramble_ref_prep/'
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
	eh_variant_catalog
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
	scramble_mode = params.scramble ? params.scramble.toString().trim().toLowerCase() : "no"
	scramble_vcf = null
	if ( scramble_mode == "yes" ) {
		scramble_ref_meta = ref_fa.map { ref_tuple -> tuple([id: 'reference'], ref_tuple[0], ref_tuple[1]) }
		scramble_ref_bundle = scramble_ref_prep(scramble_ref_meta)
		scramble_cluster_input = out_prefix.combine(bam).map { prefix, bam_tuple ->
			tuple([id: prefix], bam_tuple[0], bam_tuple[1])
		}
		scramble(scramble_cluster_input, scramble_ref_bundle.out.ref)
		scramble_vcf = scramble.out.vcf.map { meta, vcf -> vcf }
	}

	cnvnator_mode = params.cnvnator ? params.cnvnator.toString().trim().toLowerCase() : "yes"
	if ( cnvnator_mode != "no" && scramble_mode == "yes" ) {
		normalize_shortread_alignment(bam,out_prefix,ref_fa)
		CNVnator(normalize_shortread_alignment.out,out_prefix,ref_fa,params.cnvnator_bin_size)
		sv_merge_inputs = Manta.out.combine(CNVnator.out.vcf).combine(scramble_vcf).map { combined_vcfs ->
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
	else if ( scramble_mode == "yes" ) {
		sv_merge_inputs = Manta.out.combine(scramble_vcf).map { combined_vcfs ->
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
	ExpansionHunter(bam,out_prefix,ref_fa,eh_variant_catalog)
	eh_filter(out_prefix,ExpansionHunter.out)
		sv_prio(out_prefix,PhenoSV.out,ANNOVAR_SV.out,hpo,inheritance_mode,include_clinvar_report,allow_unphased_comphet)
}
