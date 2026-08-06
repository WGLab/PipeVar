
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
include { DENOVO_SNV_VCF_FILTER_CORE } from '../denovo_snv_vcf_filter_core'



// CSV batch: short-read SNP-only prioritization path with DeepVariant.
workflow INPUT_CSV_NGS_SNP {
	take:
	input_bam
	input_meta
	ref_fa
	gatk_ref_fa
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
	// Family comparison must see unrestricted cohort calls. Phenotype extraction is
	// delayed until the helper has reduced the cohort to probands.
        bam_for_proband_tasks = input_bam_with_bam
        if ( denovo_filter == "yes" ) {
                caller_regions=input_bam_with_bam.map { out_prefix, bam_file, bai_file -> tuple(out_prefix, []) }
                if ( caller_mode == "haplotypecaller" ) {
                        multi_prep_gatk_result=multi_prep_gatk(input_bam_with_bam)
                        haplotypecaller_input=multi_prep_gatk_result.join(caller_regions, failOnMismatch: true, failOnDuplicate: true)
                        snp_result=multi_haplotypecaller(haplotypecaller_input,gatk_ref_fa)
                }
                else {
                        deepvariant_input=input_bam_with_bam.join(caller_regions, failOnMismatch: true, failOnDuplicate: true)
                        snp_result=multi_deepvariant(deepvariant_input,ref_fa)
                }
                denovo_result=DENOVO_SNV_VCF_FILTER_CORE(snp_result,denovo_pedigree,denovo_role_column,denovo_family_column,denovo_vcf_sample_column,denovo_exclude_contigs)
                snp_for_annotation=denovo_result.records
                proband_keys=denovo_result.records.map { out_prefix, vcf_file -> tuple(out_prefix, true) }
                input_bam_no_bam=input_bam_no_bam.join(proband_keys, failOnDuplicate: true).map { out_prefix, note_file, marker -> tuple(out_prefix, note_file) }
                bam_for_proband_tasks=input_bam_with_bam.join(proband_keys, failOnDuplicate: true).map { out_prefix, bam_file, bai_file, marker -> tuple(out_prefix, bam_file, bai_file) }
        }
        else {
                if ( is_note == "yes" ) {
                        if ( params.phenotype_extractor.toString().trim().toLowerCase() == "phenogpt2" ) {
                                input_bam_no_bam=multi_phenogpt2(input_bam_no_bam)
                        }
                        else {
                                input_bam_no_bam=multi_phenotagger(input_bam_no_bam)
                        }
                }
                phen2gene_result=multi_phen2gene(input_bam_no_bam)
                caller_regions=input_bam_with_bam.map { out_prefix, bam_file, bai_file -> tuple(out_prefix, []) }
                if ( target == "yes" ) {
                        phen2_gene_bed=multi_phen2gene_filter(phen2gene_result,ref_fa,phen2gene_top_n)
                        caller_regions=phen2_gene_bed
                }
                if ( caller_mode == "haplotypecaller" ) {
                        multi_prep_gatk_result=multi_prep_gatk(input_bam_with_bam)
                        haplotypecaller_input=multi_prep_gatk_result.join(caller_regions, failOnMismatch: true, failOnDuplicate: true)
                        snp_result=multi_haplotypecaller(haplotypecaller_input,gatk_ref_fa)
                }
                else {
                        deepvariant_input=input_bam_with_bam.join(caller_regions, failOnMismatch: true, failOnDuplicate: true)
                        snp_result=multi_deepvariant(deepvariant_input,ref_fa)
                }
                snp_for_annotation=snp_result
        }
        if ( denovo_filter == "yes" ) {
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
                }
        }
        if ( target == "yes" ) {
                annovar_input=snp_for_annotation.join(phen2_gene_bed, failOnMismatch: true, failOnDuplicate: true).map { out_prefix, vcf_file, bed_file -> tuple(out_prefix, vcf_file, bed_file) }
        }
        else {
                annovar_input=snp_for_annotation.map { out_prefix, vcf_file -> tuple(out_prefix, vcf_file, target) }
        }
        annovar_result=multi_annovar(annovar_input)
	annovar_for_downstream=annovar_result
	annovar_result_txt=annovar_for_downstream.map { item -> tuple(item[0], item[1]) }
        join_annovar_phen2gene=annovar_result_txt.join(phen2gene_result)
        join_annovar_hpo=join_annovar_phen2gene.join(input_bam_no_bam)
        multi_eh_result=multi_expansionhunter(bam_for_proband_tasks,ref_fa,eh_variant_catalog)
        multi_eh_filter(multi_eh_result)
	rankscore_result=multi_rankscore(join_annovar_phen2gene,gnomad,rankscore_filter,rankscore_softwares,gq,ad,phen2gene_top_n)
        rankvar_result=multi_rankvar(join_annovar_hpo,gnomad,gq,ad,rankvar_filter)
        rankscore_rankvar_join=rankscore_result.join(rankvar_result)
        annovar_result_vcf=annovar_for_downstream.map { item -> tuple(item[0], item[2]) }
        snp_prio_input=rankscore_rankvar_join.join(annovar_result_vcf)
        input_bam_hpo_age=input_bam_no_bam.join(input_meta, failOnMismatch: true, failOnDuplicate: true).map { out_prefix, hpo_path, age_of_onset, sex -> tuple(out_prefix, hpo_path, age_of_onset, sex) }
        snp_prio_input_hpo=snp_prio_input.join(input_bam_hpo_age)
        multi_snp_prio(snp_prio_input_hpo,inheritance_mode,include_clinvar_report,allow_unphased_comphet)

}	
