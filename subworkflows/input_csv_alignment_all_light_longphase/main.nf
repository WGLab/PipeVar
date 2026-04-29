
include { multi_annovar } from '../../modules/multi_annovar/'
include { multi_phen2gene } from '../../modules/multi_phen2gene/'
include { multi_rankscore } from '../../modules/multi_rankscore/'
include { multi_annovar_sv } from '../../modules/multi_annovar_sv/'
include { multi_survivor } from '../../modules/multi_survivor/'
include { multi_phenosv } from '../../modules/multi_phenosv/'
include { multi_rankvar } from '../../modules/multi_rankvar/'
include { multi_sniffles } from '../../modules/multi_sniffles/'
include { multi_nanorepeat } from '../../modules/multi_nanorepeat/'
include { multi_nanocaller } from '../../modules/multi_nanocaller/'
include { multi_longphase } from '../../modules/multi_longphase/'
include { multi_phenotagger } from '../../modules/multi_phenotagger/'
include { multi_phen2gene_filter } from '../../modules/multi_reduce_region_phen2gene/'



// CSV batch: light long-read full path (SNP + SV + STR) with NanoCaller and longphase prioritization.
workflow INPUT_CSV_ALIGNMENT_ALL_LIGHT_LONGPHASE {
	take:
	input_bam
	input_age
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
                input_bam_no_bam=multi_phenotagger(input_bam_no_bam)
        }
	phen2gene_result=multi_phen2gene(input_bam_no_bam)
        if ( target == "yes" ) {
                phen2_gene_bed=multi_phen2gene_filter(phen2gene_result,ref_fa,phen2gene_top_n)
                nanocaller_result=multi_nanocaller(input_bam_with_bam,ref_fa,phen2_gene_bed)
        }
        else {
                nanocaller_result=multi_nanocaller(input_bam_with_bam,ref_fa,target)
        }
        if ( target == "yes" ) {
                annovar_input=nanocaller_result.join(phen2_gene_bed).map { out_prefix, vcf_file, bed_file -> tuple(out_prefix, vcf_file, bed_file) }
        }
        else {
                annovar_input=nanocaller_result.map { out_prefix, vcf_file -> tuple(out_prefix, vcf_file, target) }
        }
	annovar_result=multi_annovar(annovar_input)
	annovar_result_txt=annovar_result.map { item -> tuple(item[0], item[1]) }
	join_annovar_phen2gene=annovar_result_txt.join(phen2gene_result)
	join_annovar_hpo=join_annovar_phen2gene.join(input_bam_no_bam)
	rankscore_result=multi_rankscore(join_annovar_phen2gene,gnomad,rankscore_filter,rankscore_softwares,gq,phen2gene_top_n)
	rankvar_result=multi_rankvar(join_annovar_hpo,gnomad,gq,ad,rankvar_filter)
	sniffles_result=multi_sniffles(input_bam_with_bam,ref_fa)
        if ( target == "yes" ) {
	        sniffles_result_annovar=sniffles_result.join(phen2gene_result).join(phen2_gene_bed).map { out_prefix, vcf_file, phen2gene_file, bed_file -> tuple(out_prefix, vcf_file, phen2gene_file, bed_file, "called") }
        }
        else {
	        sniffles_result_annovar=sniffles_result.join(phen2gene_result).map { out_prefix, vcf_file, phen2gene_file -> tuple(out_prefix, vcf_file, phen2gene_file, target, "called") }
        }
	annovar_sv_result=multi_annovar_sv(sniffles_result_annovar)
	survivor_result=multi_survivor(annovar_sv_result)
	phenosv_input=survivor_result.join(input_bam_no_bam)
	phenosv_result=multi_phenosv(phenosv_input)
	multi_nanorepeat(input_bam_with_bam,ref_fa)
	annovar_join=annovar_result.map { item -> tuple(item[0], item[2]) }
        annovar_sv_join=annovar_sv_result.map { item -> tuple(item[0], item[1]) }
	join_vcf_bam=annovar_join.join(input_bam_with_bam)
	join_vcf_bam_sv=annovar_sv_join.join(join_vcf_bam)
	join_vcf_bam_phenosv=phenosv_result.join(join_vcf_bam_sv)
	join_vcf_bam_rankscore=rankscore_result.join(join_vcf_bam_phenosv)
	join_vcf_bam_rankvar=rankvar_result.join(join_vcf_bam_rankscore)
	input_bam_hpo_age=input_bam_no_bam.join(input_age).map { out_prefix, hpo_path, age_of_onset -> tuple(out_prefix, hpo_path, age_of_onset) }
	join_vcf_bam_rankvar_hpo=join_vcf_bam_rankvar.join(input_bam_hpo_age)
	multi_longphase(join_vcf_bam_rankvar_hpo,ref_fa,inheritance_mode,include_clinvar_report,allow_unphased_comphet)

}	
