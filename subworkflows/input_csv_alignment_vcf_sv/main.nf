
include { multi_annovar } from '../../modules/multi_annovar/'
include { multi_phen2gene } from '../../modules/multi_phen2gene/'
include { multi_rankscore } from '../../modules/multi_rankscore/'
include { multi_annovar_sv } from '../../modules/multi_annovar_sv/'
include { multi_survivor } from '../../modules/multi_survivor/'
include { multi_phenosv } from '../../modules/multi_phenosv/'
include { multi_rankvar } from '../../modules/multi_rankvar/'
include { multi_sniffles } from '../../modules/multi_sniffles/'
include { multi_manta } from '../../modules/multi_manta/'
include { multi_clair3 } from '../../modules/multi_clair3/'
include { multi_nanorepeat } from '../../modules/multi_nanorepeat/'
include { multi_expansionhunter } from '../../modules/multi_expansionhunter/'
include { multi_eh_filter } from '../../modules/multi_eh_filter/'
include { multi_deepvariant } from '../../modules/multi_deepvariant/'
include { multi_nanocaller } from '../../modules/multi_nanocaller/'
include { multi_haplotypecaller } from '../../modules/multi_haplotypecaller/'
include { multi_prep_gatk } from '../../modules/multi_prep_gatk/'
include { multi_longphase } from '../../modules/multi_longphase/'
include { multi_phenotagger } from '../../modules/multi_phenotagger/'
include { multi_sv_prio } from '../../modules/multi_sv_prio/'



workflow INPUT_CSV_ALIGNMENT_VCF_SV {
	take:
	input_vcf
	ref_fa
	is_note

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
	multi_sv_prio(sv_prio_input)
}	



