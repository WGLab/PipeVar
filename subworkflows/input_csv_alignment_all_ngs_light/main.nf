
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
include { multi_haplotypecaller } from '../../modules/multi_haplotypecaller/'
include { multi_prep_gatk } from '../../modules/multi_prep_gatk/'
include { multi_phenotagger } from '../../modules/multi_phenotagger/'
include { multi_phen2gene_filter } from '../../modules/multi_reduce_region_phen2gene/'
include { multi_ngs_prio } from '../../modules/multi_ngs_prio/'
include { multi_variant_html_report } from '../../modules/variant_html_report/'


// CSV batch: short-read full path using light SNP caller (HaplotypeCaller) plus SV path.
workflow INPUT_CSV_ALIGNMENT_ALL_NGS_LIGHT {
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
	ref_fa_no_dict = ref_fa.map { ref_fa, fai_file, dict_file -> return tuple ( ref_fa, fai_file ) }
        if ( is_note == "yes" ) {
                input_bam_no_bam=multi_phenotagger(input_bam_no_bam)
        }	
        multi_prep_gatk_result=multi_prep_gatk(input_bam)
	phen2gene_result=multi_phen2gene(input_bam_no_bam)
        if ( target == "yes" ) {
                phen2_gene_bed=multi_phen2gene_filter(phen2gene_result,ref_fa,phen2gene_top_n)
                haplotypecaller_result=multi_haplotypecaller(multi_prep_gatk_result,ref_fa,phen2_gene_bed)
        }
        else {
                haplotypecaller_result=multi_haplotypecaller(multi_prep_gatk_result,ref_fa,target)
        }
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
        rankscore_result=multi_rankscore(join_annovar_phen2gene,gnomad,rankscore_filter,rankscore_softwares,gq,phen2gene_top_n)
        rankvar_result=multi_rankvar(join_annovar_hpo,gnomad,gq,ad,rankvar_filter)
        multi_eh_result=multi_expansionhunter(input_bam_with_bam,ref_fa_no_dict)
        multi_eh_filter(multi_eh_result)
        manta_result=multi_manta(input_bam_with_bam,ref_fa_no_dict)
        if ( target == "yes" ) {
	        manta_result_annovar=manta_result.join(phen2gene_result).join(phen2_gene_bed).map { out_prefix, vcf_file, phen2gene_file, bed_file -> tuple(out_prefix, vcf_file, phen2gene_file, bed_file) }
        }
        else {
	        manta_result_annovar=manta_result.join(phen2gene_result).map { out_prefix, vcf_file, phen2gene_file -> tuple(out_prefix, vcf_file, phen2gene_file, target) }
        }
	annovar_sv_result=multi_annovar_sv(manta_result_annovar)
        survivor_result=multi_survivor(annovar_sv_result)
        phenosv_input=survivor_result.join(input_bam_no_bam)
        phenosv_result=multi_phenosv(phenosv_input)
        annovar_result_vcf=annovar_result.map { item -> tuple(item[0], item[2]) }
        phenosv_annovar_snv=phenosv_result.join(annovar_result_vcf)
        sv_join=phenosv_annovar_snv.join(annovar_sv_result)
        rankscore_join=sv_join.join(rankscore_result)
	        rankvar_join=rankscore_join.join(rankvar_result)
	        input_bam_hpo_age=input_bam_no_bam.join(input_age).map { out_prefix, hpo_path, age_of_onset -> tuple(out_prefix, hpo_path, age_of_onset) }
	        rankvar_join_hpo=rankvar_join.join(input_bam_hpo_age)
	        rankvar_join_hpo_ordered=rankvar_join_hpo.map { out_prefix, sv_pathogenic, snv_vcf_path, sv_vcf_path, snv_rankscore, snv_pathogenic, snv_rankvar, hpo_path, age_of_onset ->
	            tuple(out_prefix, snv_rankvar, snv_rankscore, snv_pathogenic, sv_pathogenic, sv_vcf_path, snv_vcf_path, hpo_path, age_of_onset)
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
	multi_variant_html_report(prio_report_input)

}	
