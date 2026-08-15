
include { multi_annovar_sv } from '../../modules/multi_annovar_sv/'
include { multi_common_sv_filter } from '../../modules/multi_common_sv_filter/'
include { multi_survivor } from '../../modules/multi_survivor/'
include { multi_phenosv } from '../../modules/multi_phenosv/'
include { multi_phenotagger } from '../../modules/multi_phenotagger/'
include { multi_phenogpt2 } from '../../modules/multi_phenogpt2/'
include { multi_sv_prio } from '../../modules/multi_sv_prio/'
include { DENOVO_SV_FILTER_CORE } from '../denovo_sv_filter_core'



// CSV batch: VCF SV re-annotation and SV prioritization path.
workflow INPUT_CSV_ALIGNMENT_VCF_SV {
	take:
	input_vcf
	input_meta
	ref_fa
	is_note
	inheritance_mode
	include_clinvar_report
	allow_unphased_comphet
	denovo_filter
	denovo_pedigree
	denovo_role_column
	denovo_family_column
	denovo_vcf_sample_column
	denovo_exclude_contigs
	denovo_sv_min_reciprocal_overlap

	main:

        input_vcf_no_vcf =  input_vcf.map { out_prefix, vcf_file, note_file -> return tuple ( out_prefix,note_file ) }
        sv_for_annotation = input_vcf.map { out_prefix, vcf_file ,note_file -> tuple(out_prefix, vcf_file) }
        if ( denovo_filter == "yes" ) {
                denovo_result = DENOVO_SV_FILTER_CORE(sv_for_annotation, denovo_pedigree, denovo_role_column, denovo_family_column, denovo_vcf_sample_column, denovo_exclude_contigs, denovo_sv_min_reciprocal_overlap)
                sv_for_annotation = denovo_result.records
                proband_keys = denovo_result.records.map { out_prefix, vcf_file -> tuple(out_prefix, true) }
                input_vcf_no_vcf = input_vcf_no_vcf.join(proband_keys, failOnDuplicate: true).map { out_prefix, note_file, marker -> tuple(out_prefix, note_file) }
        }
        if ( is_note == "yes" ) {
                if ( params.phenotype_extractor.toString().trim().toLowerCase() == "phenogpt2" ) {

                        input_vcf_no_vcf=multi_phenogpt2(input_vcf_no_vcf)

                }

                else {

                        input_vcf_no_vcf=multi_phenotagger(input_vcf_no_vcf)

                }
        }
	sv_result_annovar=sv_for_annotation.map { out_prefix, vcf_file -> tuple(out_prefix, vcf_file, "called") }
	annovar_sv_result=multi_annovar_sv(sv_result_annovar)
	annovar_sv_for_downstream = annovar_sv_result
	if ( params.common_sv_filter.toString().trim().toLowerCase() == "yes" ) {
		multi_common_sv_filter(annovar_sv_for_downstream)
		annovar_sv_for_downstream = multi_common_sv_filter.out.filtered_vcf
	}
        survivor_result=multi_survivor(annovar_sv_for_downstream)
	phenosv_input=survivor_result.join(input_vcf_no_vcf, failOnMismatch: true, failOnDuplicate: true)
        phenosv_result=multi_phenosv(phenosv_input)
	sv_prio_input=phenosv_result.join(annovar_sv_for_downstream, failOnMismatch: true, failOnDuplicate: true)
	input_vcf_hpo_age=input_vcf_no_vcf.join(input_meta, failOnMismatch: true, failOnDuplicate: true).map { out_prefix, hpo_path, age_of_onset, sex -> tuple(out_prefix, hpo_path, age_of_onset, sex) }
	sv_prio_input_hpo=sv_prio_input.join(input_vcf_hpo_age, failOnMismatch: true, failOnDuplicate: true)
	multi_sv_prio(sv_prio_input_hpo,inheritance_mode,include_clinvar_report,allow_unphased_comphet)
}	
