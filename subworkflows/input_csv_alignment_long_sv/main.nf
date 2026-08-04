
include { multi_phen2gene } from '../../modules/multi_phen2gene/'
include { multi_annovar_sv } from '../../modules/multi_annovar_sv/'
include { multi_common_sv_filter } from '../../modules/multi_common_sv_filter/'
include { multi_survivor } from '../../modules/multi_survivor/'
include { multi_phenosv } from '../../modules/multi_phenosv/'
include { multi_sniffles } from '../../modules/multi_sniffles/'
include { multi_nanorepeat } from '../../modules/multi_nanorepeat/'
include { multi_phenotagger } from '../../modules/multi_phenotagger/'
include { multi_phenogpt2 } from '../../modules/multi_phenogpt2/'
include { multi_sv_prio } from '../../modules/multi_sv_prio/'
include { DENOVO_SV_FILTER_CORE } from '../denovo_sv_filter_core'


// CSV batch: long-read SV/STR prioritization path.
workflow INPUT_CSV_ALIGNMENT_LONG_SV {
	take:
	input_bam
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

        input_bam_no_bam =  input_bam.map { out_prefix, bam_file, bai_file, note_file -> return tuple ( out_prefix,note_file ) }
	input_bam_with_bam= input_bam.map { out_prefix, bam_file, bai_file, note_file -> return tuple (out_prefix, bam_file, bai_file) }
	sniffles_result=multi_sniffles(input_bam_with_bam,ref_fa)
	sv_result=sniffles_result
	sv_for_annotation=sv_result
	bam_for_proband_tasks=input_bam_with_bam
	if ( denovo_filter == "yes" ) {
		denovo_result=DENOVO_SV_FILTER_CORE(sv_result,denovo_pedigree,denovo_role_column,denovo_family_column,denovo_vcf_sample_column,denovo_exclude_contigs,denovo_sv_min_reciprocal_overlap)
		sv_for_annotation=denovo_result.records
		proband_keys=denovo_result.records.map { out_prefix, vcf_file -> tuple(out_prefix, true) }
		input_bam_no_bam=input_bam_no_bam.join(proband_keys, failOnDuplicate: true).map { out_prefix, note_file, marker -> tuple(out_prefix, note_file) }
		bam_for_proband_tasks=input_bam_with_bam.join(proband_keys, failOnDuplicate: true).map { out_prefix, bam_file, bai_file, marker -> tuple(out_prefix, bam_file, bai_file) }
	}
        if ( is_note == "yes" ) {
                if ( params.phenotype_extractor.toString().trim().toLowerCase() == "phenogpt2" ) {
                        input_bam_no_bam=multi_phenogpt2(input_bam_no_bam)
                }
                else {
                        input_bam_no_bam=multi_phenotagger(input_bam_no_bam)
                }
	}
	phen2gene_result=multi_phen2gene(input_bam_no_bam)
	sniffles_result_annovar=sv_for_annotation.join(phen2gene_result, failOnMismatch: true, failOnDuplicate: true).map { out_prefix, vcf_file, phen2gene_file -> tuple(out_prefix, vcf_file, phen2gene_file, "null", "called") }
	annovar_sv_result=multi_annovar_sv(sniffles_result_annovar)
	annovar_sv_for_downstream = annovar_sv_result
	if ( params.common_sv_filter.toString().trim().toLowerCase() == "yes" ) {
		multi_common_sv_filter(annovar_sv_for_downstream)
		annovar_sv_for_downstream = multi_common_sv_filter.out.filtered_vcf
	}
	survivor_result=multi_survivor(annovar_sv_for_downstream)
	phenosv_input=survivor_result.join(input_bam_no_bam)
	phenosv_result=multi_phenosv(phenosv_input)
	multi_nanorepeat(bam_for_proband_tasks,ref_fa)
        sv_prio_input=phenosv_result.join(annovar_sv_for_downstream)
        input_bam_hpo_age=input_bam_no_bam.join(input_meta, failOnMismatch: true, failOnDuplicate: true).map { out_prefix, hpo_path, age_of_onset, sex -> tuple(out_prefix, hpo_path, age_of_onset, sex) }
        sv_prio_input_hpo=sv_prio_input.join(input_bam_hpo_age)
        multi_sv_prio(sv_prio_input_hpo,inheritance_mode,include_clinvar_report,allow_unphased_comphet)

}	
