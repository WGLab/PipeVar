
include { multi_annovar } from '../../modules/multi_annovar/'
include { multi_phen2gene } from '../../modules/multi_phen2gene/'
include { multi_rankscore } from '../../modules/multi_rankscore/'
include { multi_rankvar } from '../../modules/multi_rankvar/'
include { multi_haplotypecaller } from '../../modules/multi_haplotypecaller/'
include { multi_prep_gatk } from '../../modules/multi_prep_gatk/'
include { multi_phenotagger } from '../../modules/multi_phenotagger/'
include { multi_phenogpt2 } from '../../modules/multi_phenogpt2/'
include { multi_phen2gene_filter } from '../../modules/multi_reduce_region_phen2gene/'



// CSV batch: short-read SNP-only light path with HaplotypeCaller.
workflow INPUT_CSV_ALIGNMENT_NGS_SNP_LIGHT {
	take:
	input_bam
	input_meta
	ref_fa
	rankscore_filter
	rankscore_softwares
	phen2gene_top_n
	gnomad
	gq
	ad

	rankvar_filter
	is_note
	target
	inheritance_mode
	include_clinvar_report
	allow_unphased_comphet

	main:

	
	input_bam_no_bam =  input_bam.map { out_prefix, bam_file, bai_file, note_file -> return tuple ( out_prefix,note_file ) }
        input_bam_with_bam= input_bam.map { out_prefix, bam_file, bai_file, note_file -> return tuple (out_prefix, bam_file, bai_file) }
        if ( is_note == "yes" ) {
                if ( params.phenotype_extractor.toString().trim().toLowerCase() == "phenogpt2" ) {

                        input_bam_no_bam=multi_phenogpt2(input_bam_no_bam)

                }

                else {

                        input_bam_no_bam=multi_phenotagger(input_bam_no_bam)

                }
        }
        multi_prep_gatk_result=multi_prep_gatk(input_bam_with_bam)
	phen2gene_result=multi_phen2gene(input_bam_no_bam)
        caller_regions=input_bam_with_bam.map { out_prefix, bam_file, bai_file -> tuple(out_prefix, []) }
        if ( target == "yes" ) {
                phen2_gene_bed=multi_phen2gene_filter(phen2gene_result,ref_fa,phen2gene_top_n)
                caller_regions=phen2_gene_bed
        }	
        haplotypecaller_input=multi_prep_gatk_result.join(caller_regions, failOnMismatch: true, failOnDuplicate: true)
        haplotypecaller_result=multi_haplotypecaller(haplotypecaller_input,ref_fa)
        if ( target == "yes" ) {
                annovar_input=haplotypecaller_result.join(phen2_gene_bed).map { out_prefix, vcf_file, bed_file -> tuple(out_prefix, vcf_file, bed_file) }
        }
        else {
                annovar_input=haplotypecaller_result.map { out_prefix, vcf_file -> tuple(out_prefix, vcf_file, target) }
        }
        annovar_result=multi_annovar(annovar_input)
	annovar_result_txt=annovar_result.map { item -> tuple(item[0], item[1]) }
        join_annovar_phen2gene=annovar_result_txt.join(phen2gene_result)
        join_annovar_hpo=join_annovar_phen2gene.join(input_bam_no_bam)
        rankscore_result=multi_rankscore(join_annovar_phen2gene,gnomad,rankscore_filter,rankscore_softwares,gq,ad,phen2gene_top_n)
        rankvar_result=multi_rankvar(join_annovar_hpo,gnomad,gq,ad,rankvar_filter)
        rankscore_rankvar_join=rankscore_result.join(rankvar_result)
        annovar_result_vcf=annovar_result.map { item -> tuple(item[0], item[2]) }
        snp_prio_input=rankscore_rankvar_join.join(annovar_result_vcf)
        input_bam_hpo_age=input_bam_no_bam.join(input_meta).map { out_prefix, hpo_path, age_of_onset, sex -> tuple(out_prefix, hpo_path, age_of_onset, sex) }
        snp_prio_input_hpo=snp_prio_input.join(input_bam_hpo_age)
        multi_snp_prio(snp_prio_input_hpo,inheritance_mode,include_clinvar_report,allow_unphased_comphet)

}	
