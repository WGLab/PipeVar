
include { multi_annovar } from '../../modules/multi_annovar/'
include { multi_phen2gene } from '../../modules/multi_phen2gene/'
include { multi_rankscore } from '../../modules/multi_rankscore/'
include { multi_annovar_sv } from '../../modules/multi_annovar_sv/'
include { multi_survivor } from '../../modules/multi_survivor/'
include { multi_phenosv } from '../../modules/multi_phenosv/'
include { multi_rankvar } from '../../modules/multi_rankvar/'
include { multi_manta } from '../../modules/multi_manta/'
include { multi_expansionhunter } from '../../modules/multi_expansionhunter/'
include { multi_eh_filter } from '../../modules/multi_eh_filter/'
include { multi_deepvariant } from '../../modules/multi_deepvariant/'
include { multi_haplotypecaller } from '../../modules/multi_haplotypecaller/'
include { multi_prep_gatk } from '../../modules/multi_prep_gatk/'
include { multi_phenotagger } from '../../modules/multi_phenotagger/'
include { multi_phen2gene_filter } from '../../modules/multi_reduce_region_phen2gene/'
include { multi_ngs_prio } from '../../modules/multi_ngs_prio/'


// CSV batch: short-read full path (SNP + SV + STR) with DeepVariant and Manta.
workflow INPUT_CSV_ALIGNMENT_ALL_NGS {
	take:
	input_bam
	ref_fa
	rankscore_filter
	phen2gene_top_n
	gnomad
	gq
	ad

	rankvar_filter
	is_note
	target
	caller_mode

	main:

	input_bam_no_bam =  input_bam.map { out_prefix, bam_file, bai_file, note_file -> return tuple ( out_prefix,note_file ) }
        input_bam_with_bam= input_bam.map { out_prefix, bam_file, bai_file, note_file -> return tuple (out_prefix, bam_file, bai_file) }	
	if ( is_note == "yes" ) {
		input_bam_no_bam=multi_phenotagger(input_bam_no_bam)
	}
	phen2gene_result=multi_phen2gene(input_bam_no_bam)
        if ( target == "yes" ) {
                phen2_gene_bed=multi_phen2gene_filter(phen2gene_result,ref_fa,phen2gene_top_n)
                if ( caller_mode == "haplotypecaller" ) {
                        multi_prep_gatk_result=multi_prep_gatk(input_bam_with_bam)
                        snp_result=multi_haplotypecaller(multi_prep_gatk_result,ref_fa,phen2_gene_bed)
                }
                else {
                        snp_result=multi_deepvariant(input_bam_with_bam,ref_fa,phen2_gene_bed)
                }
        }
        else {
                if ( caller_mode == "haplotypecaller" ) {
                        multi_prep_gatk_result=multi_prep_gatk(input_bam_with_bam)
                        snp_result=multi_haplotypecaller(multi_prep_gatk_result,ref_fa,target)
                }
                else {
                        snp_result=multi_deepvariant(input_bam_with_bam,ref_fa,target)
                }
        }
        annovar_result=multi_annovar(snp_result)
	annovar_result_txt=annovar_result.map { item -> tuple(item[0], item[1]) }
        join_annovar_phen2gene=annovar_result_txt.join(phen2gene_result)
        join_annovar_hpo=join_annovar_phen2gene.join(input_bam_no_bam)
        rankscore_result=multi_rankscore(join_annovar_phen2gene,gnomad,rankscore_filter,gq,phen2gene_top_n)
        rankvar_result=multi_rankvar(join_annovar_hpo,gnomad,gq,ad,rankvar_filter)
        multi_eh_result=multi_expansionhunter(input_bam_with_bam,ref_fa)
        multi_eh_filter(multi_eh_result)
        manta_result=multi_manta(input_bam_with_bam,ref_fa)
	manta_result_annovar=manta_result.join(input_bam_no_bam)
	annovar_sv_result=multi_annovar_sv(manta_result_annovar)
        survivor_result=multi_survivor(annovar_sv_result)
        phenosv_input=survivor_result.join(input_bam_no_bam)
        phenosv_result=multi_phenosv(phenosv_input)
	annovar_result_vcf=annovar_result.map { item -> tuple(item[0], item[2]) }
	phenosv_annovar_snv=phenosv_result.join(annovar_result_vcf)
	sv_join=phenosv_annovar_snv.join(annovar_sv_result)
	rankscore_join=sv_join.join(rankscore_result)
	rankvar_join=rankscore_join.join(rankvar_result)
	multi_ngs_prio(rankvar_join)

}	

