
include { multi_annovar_sv } from '../../modules/multi_annovar_sv/'
include { multi_survivor } from '../../modules/multi_survivor/'
include { multi_phenosv } from '../../modules/multi_phenosv/'
include { multi_phenotagger } from '../../modules/multi_phenotagger/'
include { multi_sv_prio } from '../../modules/multi_sv_prio/'



// CSV batch: VCF SV re-annotation and SV prioritization path.
workflow INPUT_CSV_ALIGNMENT_VCF_SV {
	take:
	input_vcf
	ref_fa
	is_note
	inheritance_mode

	main:

        input_vcf_no_vcf =  input_vcf.map { out_prefix, vcf_file, note_file -> return tuple ( out_prefix,note_file ) }
        input_vcf= input_vcf.map { out_prefix, vcf_file ,note_file -> return tuple (out_prefix, vcf_file) }
        if ( is_note == "yes" ) {
                input_vcf_no_vcf=multi_phenotagger(input_vcf_no_vcf)
        }
	sv_result_annovar=input_vcf.join(input_vcf_no_vcf)
	annovar_sv_result=multi_annovar_sv(sv_result_annovar)
        survivor_result=multi_survivor(annovar_sv_result)
        phenosv_input=survivor_result.join(input_vcf_no_vcf)
        phenosv_result=multi_phenosv(phenosv_input)
	sv_prio_input=phenosv_result.join(annovar_sv_result)
	sv_prio_input_hpo=sv_prio_input.join(input_vcf_no_vcf)
	multi_sv_prio(sv_prio_input_hpo,inheritance_mode)
}	


