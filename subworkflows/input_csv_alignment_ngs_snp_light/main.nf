
include { multi_annovar } from '../../modules/multi_annovar/'
include { multi_phen2gene } from '../../modules/multi_phen2gene/'
include { multi_rankscore } from '../../modules/multi_rankscore/'
include { multi_rankvar } from '../../modules/multi_rankvar/'
include { multi_haplotypecaller } from '../../modules/multi_haplotypecaller/'
include { multi_prep_gatk } from '../../modules/multi_prep_gatk/'
include { multi_phenotagger } from '../../modules/multi_phenotagger/'
include { multi_phen2gene_filter } from '../../modules/multi_reduce_region_phen2gene/'



// CSV batch: short-read SNP-only light path with HaplotypeCaller.
workflow INPUT_CSV_ALIGNMENT_NGS_SNP_LIGHT {
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

	main:

	
	input_bam_no_bam =  input_bam.map { out_prefix, bam_file, bai_file, note_file -> return tuple ( out_prefix,note_file ) }
        input_bam_with_bam= input_bam.map { out_prefix, bam_file, bai_file, note_file -> return tuple (out_prefix, bam_file, bai_file) }
        if ( is_note == "yes" ) {
                input_bam_no_bam=multi_phenotagger(input_bam_no_bam)
        }
        multi_prep_gatk_result=multi_prep_gatk(input_bam_with_bam)
	phen2gene_result=multi_phen2gene(input_bam_no_bam)
        if ( target == "yes" ) {
                phen2_gene_bed=multi_phen2gene_filter(phen2gene_result,ref_fa,phen2gene_top_n)
                haplotypecaller_result=multi_haplotypecaller(multi_prep_gatk_result,ref_fa,phen2_gene_bed)
        }
        else {
                haplotypecaller_result=multi_haplotypecaller(multi_prep_gatk_result,ref_fa,target)
        }	
        haplotypecaller_result_annovar=haplotypecaller_result.join(input_bam_no_bam)
        annovar_result=multi_annovar(haplotypecaller_result_annovar)
	annovar_result_txt=annovar_result.map { item -> tuple(item[0], item[1]) }
        join_annovar_phen2gene=annovar_result_txt.join(phen2gene_result)
        join_annovar_hpo=join_annovar_phen2gene.join(input_bam_no_bam)
        rankscore_result=multi_rankscore(join_annovar_phen2gene,gnomad,rankscore_filter,gq,phen2gene_top_n)
        rankvar_result=multi_rankvar(join_annovar_hpo,gnomad,gq,ad,rankvar_filter)
        rankscore_rankvar_join=rankscore_result.join(rankvar_result)
        annovar_result_vcf=annovar_result.map { item -> tuple(item[0], item[2]) }
        snp_prio_input=rankscore_rankvar_join.join(annovar_result_vcf)
        snp_prio_input_hpo=snp_prio_input.join(input_bam_no_bam)
        multi_snp_prio(snp_prio_input_hpo)

}	


