
include { multi_phen2gene } from '../../modules/multi_phen2gene/'
include { multi_annovar_sv } from '../../modules/multi_annovar_sv/'
include { multi_common_sv_filter } from '../../modules/multi_common_sv_filter/'
include { multi_survivor } from '../../modules/multi_survivor/'
include { multi_phenosv } from '../../modules/multi_phenosv/'
include { multi_manta } from '../../modules/multi_manta/'
include { multi_xtea } from '../../modules/multi_xtea/'
include { multi_normalize_shortread_alignment } from '../../modules/multi_normalize_shortread_alignment/'
include { multi_cnvnator } from '../../modules/multi_cnvnator/'
include { multi_truvari_shortread_sv_merge } from '../../modules/multi_truvari_shortread_sv_merge/'
include { multi_expansionhunter } from '../../modules/multi_expansionhunter/'
include { multi_eh_filter } from '../../modules/multi_eh_filter/'
include { multi_phenotagger } from '../../modules/multi_phenotagger/'
include { multi_sv_prio } from '../../modules/multi_sv_prio/'



// CSV batch: short-read SV/STR prioritization path with Manta.
workflow INPUT_CSV_ALIGNMENT_NGS_SV {
	take:
	input_bam
	input_meta
	ref_fa
	eh_variant_catalog
	is_note
	inheritance_mode
	include_clinvar_report
	allow_unphased_comphet

	main:

        input_bam_no_bam =  input_bam.map { out_prefix, bam_file, bai_file, note_file -> return tuple ( out_prefix,note_file ) }
        input_bam_with_bam= input_bam.map { out_prefix, bam_file, bai_file, note_file -> return tuple (out_prefix, bam_file, bai_file) }
        if ( is_note == "yes" ) {
                input_bam_no_bam=multi_phenotagger(input_bam_no_bam)
        }
        phen2gene_result=multi_phen2gene(input_bam_no_bam)
        multi_eh_result=multi_expansionhunter(input_bam_with_bam,ref_fa,eh_variant_catalog)
        multi_eh_filter(multi_eh_result)
        manta_result=multi_manta(input_bam_with_bam,ref_fa)
	xtea_mode = params.xtea ? params.xtea.toString().trim().toLowerCase() : "no"
	xtea_vcf = null
	if ( xtea_mode == "yes" ) {
		xtea_input = input_bam_with_bam.map { out_prefix, bam_file, index_file ->
			tuple([id: out_prefix], bam_file, index_file)
		}
		multi_xtea(xtea_input, ref_fa)
		xtea_vcf = multi_xtea.out.vcf.map { meta, vcf -> tuple(meta.id, vcf) }
	}

	cnvnator_mode = params.cnvnator ? params.cnvnator.toString().trim().toLowerCase() : "yes"
	if ( cnvnator_mode != "no" && xtea_mode == "yes" ) {
		normalized_bam=multi_normalize_shortread_alignment(input_bam_with_bam,ref_fa)
		multi_cnvnator(normalized_bam,ref_fa,params.cnvnator_bin_size)
		merged_sv_input=manta_result.join(multi_cnvnator.out.vcf).join(xtea_vcf).map { out_prefix, manta_vcf, cnvnator_vcf, xtea_vcf ->
			tuple(out_prefix, [manta_vcf, cnvnator_vcf, xtea_vcf])
		}
		multi_truvari_shortread_sv_merge(merged_sv_input,ref_fa)
		sv_result=multi_truvari_shortread_sv_merge.out.merged_vcf
	}
	else if ( cnvnator_mode != "no" ) {
		normalized_bam=multi_normalize_shortread_alignment(input_bam_with_bam,ref_fa)
		multi_cnvnator(normalized_bam,ref_fa,params.cnvnator_bin_size)
		merged_sv_input=manta_result.join(multi_cnvnator.out.vcf).map { out_prefix, manta_vcf, cnvnator_vcf ->
			tuple(out_prefix, [manta_vcf, cnvnator_vcf])
		}
		multi_truvari_shortread_sv_merge(merged_sv_input,ref_fa)
		sv_result=multi_truvari_shortread_sv_merge.out.merged_vcf
	}
	else if ( xtea_mode == "yes" ) {
		merged_sv_input=manta_result.join(xtea_vcf).map { out_prefix, manta_vcf, xtea_vcf ->
			tuple(out_prefix, [manta_vcf, xtea_vcf])
		}
		multi_truvari_shortread_sv_merge(merged_sv_input,ref_fa)
		sv_result=multi_truvari_shortread_sv_merge.out.merged_vcf
	}
	else {
		sv_result=manta_result
	}

	sv_result_annovar=sv_result.join(phen2gene_result).map { out_prefix, vcf_file, phen2gene_file -> tuple(out_prefix, vcf_file, phen2gene_file, "null", "called") }
	annovar_sv_result=multi_annovar_sv(sv_result_annovar)
	annovar_sv_for_downstream = annovar_sv_result
	if ( params.common_sv_filter.toString().trim().toLowerCase() == "yes" ) {
		multi_common_sv_filter(annovar_sv_result)
		annovar_sv_for_downstream = multi_common_sv_filter.out.filtered_vcf
	}
	survivor_result=multi_survivor(annovar_sv_for_downstream)
        phenosv_input=survivor_result.join(input_bam_no_bam)
        phenosv_result=multi_phenosv(phenosv_input)
        sv_prio_input=phenosv_result.join(annovar_sv_for_downstream)
        input_bam_hpo_age=input_bam_no_bam.join(input_meta).map { out_prefix, hpo_path, age_of_onset, sex -> tuple(out_prefix, hpo_path, age_of_onset, sex) }
        sv_prio_input_hpo=sv_prio_input.join(input_bam_hpo_age)
        multi_sv_prio(sv_prio_input_hpo,inheritance_mode,include_clinvar_report,allow_unphased_comphet)

}	
