
include { multi_annovar } from '../../modules/multi_annovar/'
include { multi_phen2gene } from '../../modules/multi_phen2gene/'
include { multi_rankscore } from '../../modules/multi_rankscore/'
include { multi_annovar_sv } from '../../modules/multi_annovar_sv/'
include { multi_common_sv_filter } from '../../modules/multi_common_sv_filter/'
include { multi_survivor } from '../../modules/multi_survivor/'
include { multi_phenosv } from '../../modules/multi_phenosv/'
include { multi_rankvar } from '../../modules/multi_rankvar/'
include { multi_rankvar as multi_rankvar_nanocaller } from '../../modules/multi_rankvar/'
include { multi_sniffles } from '../../modules/multi_sniffles/'
include { multi_clair3 } from '../../modules/multi_clair3/'
include { multi_nanocaller } from '../../modules/multi_nanocaller/'
include { multi_nanorepeat } from '../../modules/multi_nanorepeat/'
include { multi_longphase } from '../../modules/multi_longphase/'
include { multi_phenotagger } from '../../modules/multi_phenotagger/'
include { multi_phenogpt2 } from '../../modules/multi_phenogpt2/'
include { multi_phen2gene_filter } from '../../modules/multi_reduce_region_phen2gene/'
include { multi_variant_html_report; multi_variant_html_report_with_mito } from '../../modules/variant_html_report/'
include { DENOVO_SNV_VCF_FILTER_CORE } from '../denovo_snv_vcf_filter_core'
include { DENOVO_SV_FILTER_CORE } from '../denovo_sv_filter_core'



// CSV batch: long-read full path (SNP + SV + STR) with Clair3/Sniffles and longphase prioritization.
workflow INPUT_CSV_ALIGNMENT_ALL_LONGPHASE {
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
	caller_mode
	inheritance_mode
	include_clinvar_report
	allow_unphased_comphet
	mito_tsv
	mito_mode
	denovo_filter
	denovo_pedigree
	denovo_role_column
	denovo_family_column
	denovo_vcf_sample_column
	denovo_exclude_contigs
	denovo_sv_min_reciprocal_overlap

	main:
	
	input_bam.multiMap { out_prefix, bam_file, bai_file, note_file ->
    // Define your output channels and their structures here
    no_bam: tuple(out_prefix, note_file)
    with_bam: tuple(out_prefix, bam_file, bai_file)
}.set { split_bams_ch } // 2. Assign the multiMap object to a new variable
	// 3. Access your newly structured channels
	input_bam_no_bam = split_bams_ch.no_bam
        input_bam_with_bam = split_bams_ch.with_bam
        bam_for_proband_tasks=input_bam_with_bam
        if ( denovo_filter == "yes" ) {
                caller_regions=input_bam_with_bam.map { out_prefix, bam_file, bai_file -> tuple(out_prefix, []) }
                caller_input=input_bam_with_bam.join(caller_regions, failOnMismatch: true, failOnDuplicate: true)
                if ( caller_mode == "nanocaller" ) {
                        snp_result=multi_nanocaller(caller_input,ref_fa)
                }
                else {
                        snp_result=multi_clair3(caller_input,ref_fa)
                }
                denovo_snv_result=DENOVO_SNV_VCF_FILTER_CORE(snp_result,denovo_pedigree,denovo_role_column,denovo_family_column,denovo_vcf_sample_column,denovo_exclude_contigs)
                snp_for_annotation=denovo_snv_result.records
                proband_keys=denovo_snv_result.records.map { out_prefix, vcf_file -> tuple(out_prefix, true) }
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
                caller_input=input_bam_with_bam.join(caller_regions, failOnMismatch: true, failOnDuplicate: true)
                if ( caller_mode == "nanocaller" ) {
                        snp_result=multi_nanocaller(caller_input,ref_fa)
                }
                else {
                        snp_result=multi_clair3(caller_input,ref_fa)
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
	rankscore_result=multi_rankscore(join_annovar_phen2gene,gnomad,rankscore_filter,rankscore_softwares,gq,phen2gene_top_n)
	if ( caller_mode == "nanocaller" ) {
		rankvar_result=multi_rankvar_nanocaller(join_annovar_hpo,gnomad,gq,ad,rankvar_filter)
	}
	else {
		rankvar_result=multi_rankvar(join_annovar_hpo,gnomad,gq,ad,rankvar_filter)
	}
	sniffles_result=multi_sniffles(input_bam_with_bam,ref_fa)
	sv_result=sniffles_result
	sniffles_for_phasing=sniffles_result
	sv_for_annotation=sv_result
	if ( denovo_filter == "yes" ) {
		denovo_sv_result=DENOVO_SV_FILTER_CORE(sv_result,denovo_pedigree,denovo_role_column,denovo_family_column,denovo_vcf_sample_column,denovo_exclude_contigs,denovo_sv_min_reciprocal_overlap)
		sv_for_annotation=denovo_sv_result.records
		sniffles_for_phasing=sniffles_result.join(proband_keys, failOnDuplicate: true).map { out_prefix, sniffles_vcf, marker ->
			tuple(out_prefix, sniffles_vcf)
		}
	}
        if ( target == "yes" ) {
	        sniffles_result_annovar=sv_for_annotation.join(phen2gene_result, failOnMismatch: true, failOnDuplicate: true).join(phen2_gene_bed, failOnMismatch: true, failOnDuplicate: true).map { out_prefix, vcf_file, phen2gene_file, bed_file -> tuple(out_prefix, vcf_file, phen2gene_file, bed_file, "called") }
        }
        else {
	        sniffles_result_annovar=sv_for_annotation.join(phen2gene_result, failOnMismatch: true, failOnDuplicate: true).map { out_prefix, vcf_file, phen2gene_file -> tuple(out_prefix, vcf_file, phen2gene_file, target, "called") }
        }
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
	annovar_join=annovar_for_downstream.map { item -> tuple(item[0], item[2]) }
	join_vcf_bam=annovar_join.join(bam_for_proband_tasks, failOnMismatch: true, failOnDuplicate: true)
	join_vcf_bam_sv=sniffles_for_phasing.join(annovar_sv_for_downstream, failOnMismatch: true, failOnDuplicate: true).join(join_vcf_bam, failOnMismatch: true, failOnDuplicate: true)
	join_vcf_bam_phenosv=phenosv_result.join(join_vcf_bam_sv, failOnMismatch: true, failOnDuplicate: true)
	join_vcf_bam_rankscore=rankscore_result.join(join_vcf_bam_phenosv)
	join_vcf_bam_rankvar=rankvar_result.join(join_vcf_bam_rankscore)
	input_bam_hpo_age=input_bam_no_bam.join(input_meta, failOnMismatch: true, failOnDuplicate: true).map { out_prefix, hpo_path, age_of_onset, sex -> tuple(out_prefix, hpo_path, age_of_onset, sex) }
	join_vcf_bam_rankvar_hpo=join_vcf_bam_rankvar.join(input_bam_hpo_age)
	multi_longphase(join_vcf_bam_rankvar_hpo,ref_fa,inheritance_mode,include_clinvar_report,allow_unphased_comphet)

	prio_report_input = multi_longphase.out
		.map { prio_vcf, prio_gene_report, haplotag_bam ->
			def prefix = prio_vcf.name.replaceFirst(/\.prio\.vcf$/, "")
			tuple(prefix, prio_vcf, prio_gene_report)
		}
		.join(
			multi_nanorepeat.out.map { repeat_tsv ->
				def prefix = repeat_tsv.name.replaceFirst(/_nanorepeat_result\.tsv$/, "")
				tuple(prefix, repeat_tsv)
			}
		)
		.map { out_prefix, prio_vcf, prio_gene_report, repeat_tsv ->
			tuple(out_prefix, prio_vcf, prio_gene_report, repeat_tsv)
		}
	if ( mito_mode == "yes" ) {
		prio_report_input_with_mito = prio_report_input
			.join(
				mito_tsv.map { out_prefix, mito_report ->
					tuple(out_prefix, mito_report)
				}
			)
			.map { out_prefix, prio_vcf, prio_gene_report, repeat_tsv, mito_report ->
				tuple(out_prefix, prio_vcf, prio_gene_report, repeat_tsv, mito_report)
			}
		multi_variant_html_report_with_mito(prio_report_input_with_mito)
	}
	else {
		multi_variant_html_report(prio_report_input)
	}

}	
