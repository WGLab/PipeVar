
include { multi_annovar } from '../../modules/multi_annovar/'
include { multi_phen2gene } from '../../modules/multi_phen2gene/'
include { multi_rankscore } from '../../modules/multi_rankscore/'
include { multi_annovar_sv } from '../../modules/multi_annovar_sv/'
include { multi_common_sv_filter } from '../../modules/multi_common_sv_filter/'
include { multi_survivor } from '../../modules/multi_survivor/'
include { multi_phenosv } from '../../modules/multi_phenosv/'
include { multi_rankvar } from '../../modules/multi_rankvar/'
include { multi_sniffles } from '../../modules/multi_sniffles/'
include { multi_cnvpytor } from '../../modules/multi_cnvpytor/'
include { multi_merge_longread_sv_callers } from '../../modules/multi_merge_longread_sv_callers/'
include { multi_clair3 } from '../../modules/multi_clair3/'
include { multi_nanocaller } from '../../modules/multi_nanocaller/'
include { multi_nanorepeat } from '../../modules/multi_nanorepeat/'
include { multi_longphase } from '../../modules/multi_longphase/'
include { multi_phenotagger } from '../../modules/multi_phenotagger/'
include { multi_phen2gene_filter } from '../../modules/multi_reduce_region_phen2gene/'
include { multi_variant_html_report; multi_variant_html_report_with_mito } from '../../modules/variant_html_report/'



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

	main:
	
	input_bam.multiMap { out_prefix, bam_file, bai_file, note_file ->
    // Define your output channels and their structures here
    no_bam: tuple(out_prefix, note_file)
    with_bam: tuple(out_prefix, bam_file, bai_file)
}.set { split_bams_ch } // 2. Assign the multiMap object to a new variable
	// 3. Access your newly structured channels
	input_bam_no_bam = split_bams_ch.no_bam
	input_bam_with_bam = split_bams_ch.with_bam
	if ( is_note == "yes" ) {
		input_bam_no_bam=multi_phenotagger(input_bam_no_bam)
	}
	phen2gene_result=multi_phen2gene(input_bam_no_bam)
        if ( target == "yes" ) {
                phen2_gene_bed=multi_phen2gene_filter(phen2gene_result,ref_fa,phen2gene_top_n)
                if ( caller_mode == "nanocaller" ) {
                        snp_result=multi_nanocaller(input_bam_with_bam,ref_fa,phen2_gene_bed)
                }
                else {
                        snp_result=multi_clair3(input_bam_with_bam,ref_fa,phen2_gene_bed)
                }
        }
        else {
                if ( caller_mode == "nanocaller" ) {
                        snp_result=multi_nanocaller(input_bam_with_bam,ref_fa,target)
                }
                else {
                        snp_result=multi_clair3(input_bam_with_bam,ref_fa,target)
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
	rankscore_result=multi_rankscore(join_annovar_phen2gene,gnomad,rankscore_filter,rankscore_softwares,gq,phen2gene_top_n)
	rankvar_result=multi_rankvar(join_annovar_hpo,gnomad,gq,ad,rankvar_filter)
	sniffles_result=multi_sniffles(input_bam_with_bam,ref_fa)
	cnvpytor_mode = params.cnvpytor ? params.cnvpytor.toString().trim().toLowerCase() : "no"
	cnvpytor_baf_mode = params.cnvpytor_baf ? params.cnvpytor_baf.toString().trim().toLowerCase() : "yes"
	if ( cnvpytor_mode == "yes" ) {
		cnvpytor_input=input_bam_with_bam.join(snp_result).map { out_prefix, bam_file, bai_file, snp_vcf ->
			tuple(out_prefix, bam_file, bai_file, snp_vcf)
		}
		cnvpytor_result=multi_cnvpytor(cnvpytor_input,ref_fa,Channel.value(cnvpytor_baf_mode),Channel.value(params.cnvpytor_bin_sizes),Channel.value(params.cnvpytor_primary_bin),Channel.value(params.cnvpytor_min_size))
		merged_sv_input=sniffles_result.join(cnvpytor_result.vcf).map { out_prefix, sniffles_vcf, cnvpytor_vcf ->
			tuple(out_prefix, [sniffles_vcf, cnvpytor_vcf])
		}
		sv_result=multi_merge_longread_sv_callers(merged_sv_input)
	}
	else {
		sv_result=sniffles_result
	}
        if ( target == "yes" ) {
	        sniffles_result_annovar=sv_result.join(phen2gene_result).join(phen2_gene_bed).map { out_prefix, vcf_file, phen2gene_file, bed_file -> tuple(out_prefix, vcf_file, phen2gene_file, bed_file, "called") }
        }
        else {
	        sniffles_result_annovar=sv_result.join(phen2gene_result).map { out_prefix, vcf_file, phen2gene_file -> tuple(out_prefix, vcf_file, phen2gene_file, target, "called") }
        }
	annovar_sv_result=multi_annovar_sv(sniffles_result_annovar)
	annovar_sv_for_downstream = annovar_sv_result
	if ( params.common_sv_filter.toString().trim().toLowerCase() == "yes" ) {
		multi_common_sv_filter(annovar_sv_result)
		annovar_sv_for_downstream = multi_common_sv_filter.out.filtered_vcf
	}
	survivor_result=multi_survivor(annovar_sv_for_downstream)
	phenosv_input=survivor_result.join(input_bam_no_bam)
	phenosv_result=multi_phenosv(phenosv_input)
	multi_nanorepeat(input_bam_with_bam,ref_fa)
	annovar_join=annovar_result.map { item -> tuple(item[0], item[2]) }
	annovar_sv_join=annovar_sv_for_downstream.map { item -> tuple(item[0], item[1]) }
	join_vcf_bam=annovar_join.join(input_bam_with_bam)
	join_vcf_bam_sv=annovar_sv_join.join(join_vcf_bam)
	join_vcf_bam_phenosv=phenosv_result.join(join_vcf_bam_sv)
	join_vcf_bam_rankscore=rankscore_result.join(join_vcf_bam_phenosv)
	join_vcf_bam_rankvar=rankvar_result.join(join_vcf_bam_rankscore)
	input_bam_hpo_age=input_bam_no_bam.join(input_meta).map { out_prefix, hpo_path, age_of_onset, sex -> tuple(out_prefix, hpo_path, age_of_onset, sex) }
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
