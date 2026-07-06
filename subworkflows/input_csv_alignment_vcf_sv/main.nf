
include { multi_annovar_sv } from '../../modules/multi_annovar_sv/'
include { multi_common_sv_filter } from '../../modules/multi_common_sv_filter/'
include { multi_survivor } from '../../modules/multi_survivor/'
include { multi_phenosv } from '../../modules/multi_phenosv/'
include { multi_phen2gene } from '../../modules/multi_phen2gene/'
include { multi_phenotagger } from '../../modules/multi_phenotagger/'
include { multi_phenogpt2 } from '../../modules/multi_phenogpt2/'
include { multi_phen2gene_filter } from '../../modules/multi_reduce_region_phen2gene/'
include { multi_sv_prio } from '../../modules/multi_sv_prio/'



// CSV batch: VCF SV re-annotation and SV prioritization path.
workflow INPUT_CSV_ALIGNMENT_VCF_SV {
	take:
	input_vcf
	input_meta
	ref_fa
	phen2gene_top_n
	is_note
	target
	inheritance_mode
	include_clinvar_report
	allow_unphased_comphet

	main:

        input_vcf_no_vcf =  input_vcf.map { out_prefix, vcf_file, note_file -> return tuple ( out_prefix,note_file ) }
        input_vcf= input_vcf.map { out_prefix, vcf_file ,note_file -> return tuple (out_prefix, vcf_file) }
        if ( is_note == "yes" ) {
                if ( params.phenotype_extractor.toString().trim().toLowerCase() == "phenogpt2" ) {

                        input_vcf_no_vcf=multi_phenogpt2(input_vcf_no_vcf)

                }

                else {

                        input_vcf_no_vcf=multi_phenotagger(input_vcf_no_vcf)

                }
        }
        phen2gene_result=multi_phen2gene(input_vcf_no_vcf)
        if ( target == "yes" ) {
	        phen2_gene_bed=multi_phen2gene_filter(phen2gene_result,ref_fa,phen2gene_top_n)
	        sv_result_annovar=input_vcf.join(phen2gene_result).join(phen2_gene_bed).map { out_prefix, vcf_file, phen2gene_file, bed_file -> tuple(out_prefix, vcf_file, phen2gene_file, bed_file, "called") }
        }
        else {
	        sv_result_annovar=input_vcf.join(phen2gene_result).map { out_prefix, vcf_file, phen2gene_file -> tuple(out_prefix, vcf_file, phen2gene_file, target, "called") }
        }
	annovar_sv_result=multi_annovar_sv(sv_result_annovar)
	annovar_sv_for_downstream = annovar_sv_result
	if ( params.common_sv_filter.toString().trim().toLowerCase() == "yes" ) {
		multi_common_sv_filter(annovar_sv_result)
		annovar_sv_for_downstream = multi_common_sv_filter.out.filtered_vcf
	}
        survivor_result=multi_survivor(annovar_sv_for_downstream)
        phenosv_input=survivor_result.join(input_vcf_no_vcf)
        phenosv_result=multi_phenosv(phenosv_input)
	sv_prio_input=phenosv_result.join(annovar_sv_for_downstream)
	input_vcf_hpo_age=input_vcf_no_vcf.join(input_meta).map { out_prefix, hpo_path, age_of_onset, sex -> tuple(out_prefix, hpo_path, age_of_onset, sex) }
	sv_prio_input_hpo=sv_prio_input.join(input_vcf_hpo_age)
	multi_sv_prio(sv_prio_input_hpo,inheritance_mode,include_clinvar_report,allow_unphased_comphet)
}	
