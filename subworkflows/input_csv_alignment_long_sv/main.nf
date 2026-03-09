
include { multi_phen2gene } from '../../modules/multi_phen2gene/'
include { multi_annovar_sv } from '../../modules/multi_annovar_sv/'
include { multi_survivor } from '../../modules/multi_survivor/'
include { multi_phenosv } from '../../modules/multi_phenosv/'
include { multi_sniffles } from '../../modules/multi_sniffles/'
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
	sniffles_result_annovar=sniffles_result.join(phen2gene_result).map { out_prefix, vcf_file, phen2gene_file -> tuple(out_prefix, vcf_file, phen2gene_file, "null") }
	annovar_sv_result=multi_annovar_sv(sniffles_result_annovar)
	survivor_result=multi_survivor(annovar_sv_result)
	phenosv_input=survivor_result.join(input_bam_no_bam)
	phenosv_result=multi_phenosv(phenosv_input)
	multi_nanorepeat(input_bam_with_bam,ref_fa)
        sv_prio_input=phenosv_result.join(annovar_sv_result)
        input_bam_hpo_age=input_bam_no_bam.join(input_age).map { out_prefix, hpo_path, age_of_onset -> tuple(out_prefix, hpo_path, age_of_onset) }
        sv_prio_input_hpo=sv_prio_input.join(input_bam_hpo_age)
        multi_sv_prio(sv_prio_input_hpo,inheritance_mode,include_clinvar_report,allow_unphased_comphet)

}	
