
include { multi_annovar } from '../../modules/multi_annovar/'
include { multi_phen2gene } from '../../modules/multi_phen2gene/'
include { multi_rankscore } from '../../modules/multi_rankscore/'
include { multi_rankvar } from '../../modules/multi_rankvar/'
include { multi_phenotagger } from '../../modules/multi_phenotagger/'
include { multi_snp_prio } from '../../modules/multi_snp_prio/'


// CSV batch: VCF SNP re-annotation and SNP prioritization path.
workflow INPUT_CSV_ALIGNMENT_VCF_SNP {
	take:
	input_vcf
	ref_fa
	rankscore_filter
	phen2gene_top_n
	gnomad
	gq
	ad

	rankvar_filter
	is_note

	main:

        input_vcf_no_vcf =  input_vcf.map { out_prefix, vcf_file, note_file -> return tuple ( out_prefix,note_file ) }
        input_vcf= input_vcf.map { out_prefix, vcf_file ,note_file -> return tuple (out_prefix, vcf_file) }
        if ( is_note == "yes" ) {
                input_vcf_no_vcf=multi_phenotagger(input_vcf_no_vcf)
        }
	
        annovar_result=multi_annovar(input_vcf)
        annovar_result_txt=annovar_result.map { item -> tuple(item[0], item[1]) }
        phen2gene_result=multi_phen2gene(input_vcf_no_vcf)
        join_annovar_phen2gene=annovar_result_txt.join(phen2gene_result)
        join_annovar_hpo=join_annovar_phen2gene.join(input_vcf_no_vcf)
        rankscore_result=multi_rankscore(join_annovar_phen2gene,gnomad,rankscore_filter,gq,phen2gene_top_n)
        rankvar_result=multi_rankvar(join_annovar_hpo,gnomad,gq,ad,rankvar_filter)
	rankscore_rankvar_join=rankvar_result.join(rankscore_result)
	annovar_result_vcf=annovar_result.map { item -> tuple(item[0], item[2]) }
	snp_prio_input=rankscore_rankvar_join.join(annovar_result_vcf)
	snp_prio_input_hpo=snp_prio_input.join(input_vcf_no_vcf)
	multi_snp_prio(snp_prio_input_hpo)
}	


