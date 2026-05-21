

include { SURVIVOR } from '../../modules/survivor/'
include { PhenoSV } from '../../modules/phenosv/'
include { Phen2gene } from '../../modules/phen2gene/'
include { Manta } from '../../modules/manta/'
include { xtea } from '../../modules/xtea/'
include { normalize_shortread_alignment } from '../../modules/normalize_shortread_alignment/'
include { CNVnator } from '../../modules/cnvnator/'
include { truvari_shortread_sv_merge } from '../../modules/truvari_shortread_sv_merge/'
include { ANNOVAR_SV } from '../../modules/annovar_sv/'
include { common_sv_filter } from '../../modules/common_sv_filter/'
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
	xtea_mode = params.xtea ? params.xtea.toString().trim().toLowerCase() : "no"
	xtea_vcf = null
	if ( xtea_mode == "yes" ) {
		xtea_input = out_prefix.combine(bam).map { prefix, bam_tuple ->
			tuple([id: prefix], bam_tuple[0], bam_tuple[1])
		}
		xtea(xtea_input, ref_fa)
		xtea_vcf = xtea.out.vcf.map { meta, vcf -> vcf }
	}

	cnvnator_mode = params.cnvnator ? params.cnvnator.toString().trim().toLowerCase() : "yes"
	if ( cnvnator_mode != "no" && xtea_mode == "yes" ) {
		normalize_shortread_alignment(bam,out_prefix,ref_fa)
		CNVnator(normalize_shortread_alignment.out,out_prefix,ref_fa,params.cnvnator_bin_size)
		sv_merge_inputs = Manta.out.combine(CNVnator.out.vcf).combine(xtea_vcf).map { combined_vcfs ->
			combined_vcfs.flatten()
		}
		truvari_shortread_sv_merge(sv_merge_inputs,ref_fa,out_prefix)
		sv_vcf=truvari_shortread_sv_merge.out.merged_vcf
	}
	else if ( cnvnator_mode != "no" ) {
		normalize_shortread_alignment(bam,out_prefix,ref_fa)
		CNVnator(normalize_shortread_alignment.out,out_prefix,ref_fa,params.cnvnator_bin_size)
		sv_merge_inputs = Manta.out.combine(CNVnator.out.vcf).map { combined_vcfs ->
			combined_vcfs.flatten()
		}
		truvari_shortread_sv_merge(sv_merge_inputs,ref_fa,out_prefix)
		sv_vcf=truvari_shortread_sv_merge.out.merged_vcf
	}
	else if ( xtea_mode == "yes" ) {
		sv_merge_inputs = Manta.out.combine(xtea_vcf).map { combined_vcfs ->
			combined_vcfs.flatten()
		}
		truvari_shortread_sv_merge(sv_merge_inputs,ref_fa,out_prefix)
		sv_vcf=truvari_shortread_sv_merge.out.merged_vcf
	}
	else {
		sv_vcf=Manta.out
	}

	ANNOVAR_SV(sv_vcf,out_prefix,Phen2gene.out,"null","called")
	annovar_sv_for_downstream = ANNOVAR_SV.out
	if ( params.common_sv_filter.toString().trim().toLowerCase() == "yes" ) {
		common_sv_filter(ANNOVAR_SV.out,out_prefix)
		annovar_sv_for_downstream = common_sv_filter.out.filtered_vcf
	}
	SURVIVOR(annovar_sv_for_downstream,out_prefix)
	PhenoSV(SURVIVOR.out,out_prefix,hpo)
	ExpansionHunter(bam,out_prefix,ref_fa,eh_variant_catalog)
	eh_filter(out_prefix,ExpansionHunter.out)
		sv_prio(out_prefix,PhenoSV.out,annovar_sv_for_downstream,hpo,inheritance_mode,include_clinvar_report,allow_unphased_comphet)
}
