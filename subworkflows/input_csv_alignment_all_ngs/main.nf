
include { multi_annovar } from '../../modules/multi_annovar/'
include { multi_phen2gene } from '../../modules/multi_phen2gene/'
include { multi_rankscore } from '../../modules/multi_rankscore/'
include { multi_annovar_sv } from '../../modules/multi_annovar_sv/'
include { multi_common_sv_filter } from '../../modules/multi_common_sv_filter/'
include { multi_survivor } from '../../modules/multi_survivor/'
include { multi_phenosv } from '../../modules/multi_phenosv/'
include { multi_rankvar } from '../../modules/multi_rankvar/'
include { multi_manta } from '../../modules/multi_manta/'
include { multi_xtea } from '../../modules/multi_xtea/'
include { multi_normalize_shortread_alignment } from '../../modules/multi_normalize_shortread_alignment/'
include { multi_cnvnator } from '../../modules/multi_cnvnator/'
include { multi_truvari_shortread_sv_merge } from '../../modules/multi_truvari_shortread_sv_merge/'
include { multi_expansionhunter } from '../../modules/multi_expansionhunter/'
include { multi_eh_filter } from '../../modules/multi_eh_filter/'
include { multi_deepvariant } from '../../modules/multi_deepvariant/'
include { multi_haplotypecaller } from '../../modules/multi_haplotypecaller/'
include { multi_prep_gatk } from '../../modules/multi_prep_gatk/'
include { multi_phenotagger } from '../../modules/multi_phenotagger/'
include { multi_phenogpt2 } from '../../modules/multi_phenogpt2/'
include { multi_phen2gene_filter } from '../../modules/multi_reduce_region_phen2gene/'
include { multi_ngs_prio } from '../../modules/multi_ngs_prio/'
include { multi_variant_html_report; multi_variant_html_report_with_mito } from '../../modules/variant_html_report/'
include { DENOVO_SNV_VCF_FILTER_CORE } from '../denovo_snv_vcf_filter_core'
include { DENOVO_SV_FILTER_CORE } from '../denovo_sv_filter_core'


// CSV batch: short-read full path (SNP + SV + STR) with DeepVariant and Manta.
workflow INPUT_CSV_ALIGNMENT_ALL_NGS {
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
	input_bam_no_bam =  input_bam.map { out_prefix, bam_file, bai_file, note_file -> return tuple ( out_prefix,note_file ) }
        input_bam_with_bam= input_bam.map { out_prefix, bam_file, bai_file, note_file -> return tuple (out_prefix, bam_file, bai_file) }	
        bam_for_proband_tasks=input_bam_with_bam
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
        rankscore_result=multi_rankscore(join_annovar_phen2gene,gnomad,rankscore_filter,rankscore_softwares,gq,phen2gene_top_n)
        rankvar_result=multi_rankvar(join_annovar_hpo,gnomad,gq,ad,rankvar_filter)
        multi_eh_result=multi_expansionhunter(bam_for_proband_tasks,ref_fa,eh_variant_catalog)
        multi_eh_filter(multi_eh_result)
        manta_result=multi_manta(input_bam_with_bam,ref_fa)
	xtea_mode = params.xtea ? params.xtea.toString().trim().toLowerCase() : "no"
	xtea_vcf = null
	if ( xtea_mode == "yes" ) {
		xtea_input = input_bam_with_bam.map { out_prefix, bam_file, index_file ->
			tuple([id: out_prefix], bam_file, index_file)
		}
		multi_xtea(xtea_input, ref_fa)
		xtea_vcf = multi_xtea.out.vcf.map { meta, vcf -> tuple(meta.id, vcf) }
	}

	cnvnator_mode = params.cnvnator ? params.cnvnator.toString().trim().toLowerCase() : "yes"
	if ( cnvnator_mode != "no" && xtea_mode == "yes" ) {
		normalized_bam=multi_normalize_shortread_alignment(input_bam_with_bam,ref_fa)
		multi_cnvnator(normalized_bam,ref_fa,params.cnvnator_bin_size)
		merged_sv_input=manta_result.join(multi_cnvnator.out.vcf).join(xtea_vcf).map { out_prefix, manta_vcf, cnvnator_vcf, xtea_vcf ->
			tuple(out_prefix, [manta_vcf, cnvnator_vcf, xtea_vcf])
		}
		multi_truvari_shortread_sv_merge(merged_sv_input,ref_fa)
		sv_result=multi_truvari_shortread_sv_merge.out.merged_vcf
	}
	else if ( cnvnator_mode != "no" ) {
		normalized_bam=multi_normalize_shortread_alignment(input_bam_with_bam,ref_fa)
		multi_cnvnator(normalized_bam,ref_fa,params.cnvnator_bin_size)
		merged_sv_input=manta_result.join(multi_cnvnator.out.vcf).map { out_prefix, manta_vcf, cnvnator_vcf ->
			tuple(out_prefix, [manta_vcf, cnvnator_vcf])
		}
		multi_truvari_shortread_sv_merge(merged_sv_input,ref_fa)
		sv_result=multi_truvari_shortread_sv_merge.out.merged_vcf
	}
	else if ( xtea_mode == "yes" ) {
		merged_sv_input=manta_result.join(xtea_vcf).map { out_prefix, manta_vcf, xtea_vcf ->
			tuple(out_prefix, [manta_vcf, xtea_vcf])
		}
		multi_truvari_shortread_sv_merge(merged_sv_input,ref_fa)
		sv_result=multi_truvari_shortread_sv_merge.out.merged_vcf
	}
	else {
		sv_result=manta_result
	}
	sv_for_annotation=sv_result
	if ( denovo_filter == "yes" ) {
		denovo_sv_result=DENOVO_SV_FILTER_CORE(sv_result,denovo_pedigree,denovo_role_column,denovo_family_column,denovo_vcf_sample_column,denovo_exclude_contigs,denovo_sv_min_reciprocal_overlap)
		sv_for_annotation=denovo_sv_result.records
	}

        if ( target == "yes" ) {
	        sv_result_annovar=sv_for_annotation.join(phen2gene_result, failOnMismatch: true, failOnDuplicate: true).join(phen2_gene_bed, failOnMismatch: true, failOnDuplicate: true).map { out_prefix, vcf_file, phen2gene_file, bed_file -> tuple(out_prefix, vcf_file, phen2gene_file, bed_file, "called") }
        }
        else {
	        sv_result_annovar=sv_for_annotation.join(phen2gene_result, failOnMismatch: true, failOnDuplicate: true).map { out_prefix, vcf_file, phen2gene_file -> tuple(out_prefix, vcf_file, phen2gene_file, target, "called") }
        }
	annovar_sv_result=multi_annovar_sv(sv_result_annovar)
	annovar_sv_for_downstream = annovar_sv_result
	if ( params.common_sv_filter.toString().trim().toLowerCase() == "yes" ) {
		multi_common_sv_filter(annovar_sv_for_downstream)
		annovar_sv_for_downstream = multi_common_sv_filter.out.filtered_vcf
	}
        survivor_result=multi_survivor(annovar_sv_for_downstream)
        phenosv_input=survivor_result.join(input_bam_no_bam)
        phenosv_result=multi_phenosv(phenosv_input)
	annovar_result_vcf=annovar_for_downstream.map { item -> tuple(item[0], item[2]) }
	phenosv_annovar_snv=phenosv_result.join(annovar_result_vcf)
	sv_join=phenosv_annovar_snv.join(annovar_sv_for_downstream)
	rankscore_join=sv_join.join(rankscore_result)
	rankvar_join=rankscore_join.join(rankvar_result)
	input_bam_hpo_age=input_bam_no_bam.join(input_meta, failOnMismatch: true, failOnDuplicate: true).map { out_prefix, hpo_path, age_of_onset, sex -> tuple(out_prefix, hpo_path, age_of_onset, sex) }
	rankvar_join_hpo=rankvar_join.join(input_bam_hpo_age)
	rankvar_join_hpo_ordered=rankvar_join_hpo.map { out_prefix, sv_pathogenic, snv_vcf_path, sv_vcf_path, snv_rankscore, snv_pathogenic, snv_rankvar, hpo_path, age_of_onset, sex ->
	    tuple(out_prefix, snv_rankvar, snv_rankscore, snv_pathogenic, sv_pathogenic, sv_vcf_path, snv_vcf_path, hpo_path, age_of_onset, sex)
	}
	multi_ngs_prio(rankvar_join_hpo_ordered,inheritance_mode,include_clinvar_report,allow_unphased_comphet)

	prio_report_input = multi_ngs_prio.out[0]
		.map { prio_vcf ->
			def prefix = prio_vcf.name.replaceFirst(/\.prio\.vcf$/, "")
			tuple(prefix, prio_vcf)
		}
		.join(
			multi_ngs_prio.out[1].map { prio_gene_report ->
				def prefix = prio_gene_report.name.replaceFirst(/\.prio_gene\.vcf$/, "")
				tuple(prefix, prio_gene_report)
			}
		)
		.join(
			multi_eh_filter.out.map { repeat_tsv ->
				def prefix = repeat_tsv.name.replaceFirst(/\.eh\.tsv$/, "")
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
