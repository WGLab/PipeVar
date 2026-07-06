
include { multi_annovar } from '../../modules/multi_annovar/'
include { multi_phen2gene } from '../../modules/multi_phen2gene/'
include { multi_rankscore } from '../../modules/multi_rankscore/'
include { multi_rankvar } from '../../modules/multi_rankvar/'
include { multi_expansionhunter } from '../../modules/multi_expansionhunter/'
include { multi_eh_filter } from '../../modules/multi_eh_filter/'
include { multi_deepvariant } from '../../modules/multi_deepvariant/'
include { multi_haplotypecaller } from '../../modules/multi_haplotypecaller/'
include { multi_prep_gatk } from '../../modules/multi_prep_gatk/'
include { multi_phenotagger } from '../../modules/multi_phenotagger/'
include { multi_phenogpt2 } from '../../modules/multi_phenogpt2/'
include { multi_phen2gene_filter } from '../../modules/multi_reduce_region_phen2gene/'
include { multi_snp_prio } from '../../modules/multi_snp_prio/'



// CSV batch: short-read SNP-only prioritization path with DeepVariant.
workflow INPUT_CSV_NGS_SNP {
	take:
	input_bam
	input_meta
	ref_fa
	eh_variant_catalog
	rankscore_filter
	rankscore_softwares
	phen2gene_top_n
	gnomad
	gq
	ad

	rankvar_filter
	is_note
	target
	caller_mode
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
        if ( target == "yes" ) {
                annovar_input=snp_result.join(phen2_gene_bed).map { out_prefix, vcf_file, bed_file -> tuple(out_prefix, vcf_file, bed_file) }
        }
        else {
                annovar_input=snp_result.map { out_prefix, vcf_file -> tuple(out_prefix, vcf_file, target) }
        }
        annovar_result=multi_annovar(annovar_input)
	annovar_result_txt=annovar_result.map { item -> tuple(item[0], item[1]) }
        join_annovar_phen2gene=annovar_result_txt.join(phen2gene_result)
        join_annovar_hpo=join_annovar_phen2gene.join(input_bam_no_bam)
        multi_eh_result=multi_expansionhunter(input_bam_with_bam,ref_fa,eh_variant_catalog)
        multi_eh_filter(multi_eh_result)
	rankscore_result=multi_rankscore(join_annovar_phen2gene,gnomad,rankscore_filter,rankscore_softwares,gq,phen2gene_top_n)
        rankvar_result=multi_rankvar(join_annovar_hpo,gnomad,gq,ad,rankvar_filter)
        rankscore_rankvar_join=rankscore_result.join(rankvar_result)
        annovar_result_vcf=annovar_result.map { item -> tuple(item[0], item[2]) }
        snp_prio_input=rankscore_rankvar_join.join(annovar_result_vcf)
        input_bam_hpo_age=input_bam_no_bam.join(input_meta).map { out_prefix, hpo_path, age_of_onset, sex -> tuple(out_prefix, hpo_path, age_of_onset, sex) }
        snp_prio_input_hpo=snp_prio_input.join(input_bam_hpo_age)
        multi_snp_prio(snp_prio_input_hpo,inheritance_mode,include_clinvar_report,allow_unphased_comphet)

}	
