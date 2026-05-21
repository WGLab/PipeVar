
include { multi_phen2gene } from '../../modules/multi_phen2gene/'
include { multi_annovar_sv } from '../../modules/multi_annovar_sv/'
include { multi_common_sv_filter } from '../../modules/multi_common_sv_filter/'
include { multi_survivor } from '../../modules/multi_survivor/'
include { multi_phenosv } from '../../modules/multi_phenosv/'
include { multi_sniffles } from '../../modules/multi_sniffles/'
include { multi_cnvpytor } from '../../modules/multi_cnvpytor/'
include { multi_merge_longread_sv_callers } from '../../modules/multi_merge_longread_sv_callers/'
include { multi_nanorepeat } from '../../modules/multi_nanorepeat/'
include { multi_phenotagger } from '../../modules/multi_phenotagger/'
include { multi_sv_prio } from '../../modules/multi_sv_prio/'


// CSV batch: long-read SV/STR prioritization path.
workflow INPUT_CSV_ALIGNMENT_LONG_SV {
	take:
	input_bam
	input_age
	ref_fa
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
	sniffles_result=multi_sniffles(input_bam_with_bam,ref_fa)
	cnvpytor_mode = params.cnvpytor ? params.cnvpytor.toString().trim().toLowerCase() : "no"
	if ( cnvpytor_mode == "yes" ) {
		cnvpytor_input = input_bam_with_bam.map { out_prefix, bam_file, bai_file -> tuple(out_prefix, bam_file, bai_file, bam_file) }
		cnvpytor_result=multi_cnvpytor(cnvpytor_input,Channel.value("no"),Channel.value(params.cnvpytor_bin_sizes),Channel.value(params.cnvpytor_primary_bin),Channel.value(params.cnvpytor_min_size))
		merged_sv_input=sniffles_result.join(cnvpytor_result.vcf).map { out_prefix, sniffles_vcf, cnvpytor_vcf ->
			tuple(out_prefix, [sniffles_vcf, cnvpytor_vcf])
		}
		sv_result=multi_merge_longread_sv_callers(merged_sv_input)
	}
	else {
		sv_result=sniffles_result
	}
	sniffles_result_annovar=sv_result.join(phen2gene_result).map { out_prefix, vcf_file, phen2gene_file -> tuple(out_prefix, vcf_file, phen2gene_file, "null", "called") }
	annovar_sv_result=multi_annovar_sv(sniffles_result_annovar)
	annovar_sv_for_downstream = annovar_sv_result
	if ( params.common_sv_filter.toString().trim().toLowerCase() == "yes" ) {
		multi_common_sv_filter(annovar_sv_result)
		annovar_sv_for_downstream = multi_common_sv_filter.out.filtered_vcf
	}
	survivor_result=multi_survivor(annovar_sv_for_downstream)
	phenosv_input=survivor_result.join(input_bam_no_bam)
	phenosv_result=multi_phenosv(phenosv_input)
	multi_nanorepeat(input_bam_with_bam,ref_fa)
        sv_prio_input=phenosv_result.join(annovar_sv_for_downstream)
        input_bam_hpo_age=input_bam_no_bam.join(input_age).map { out_prefix, hpo_path, age_of_onset -> tuple(out_prefix, hpo_path, age_of_onset) }
        sv_prio_input_hpo=sv_prio_input.join(input_bam_hpo_age)
        multi_sv_prio(sv_prio_input_hpo,inheritance_mode,include_clinvar_report,allow_unphased_comphet)

}	
