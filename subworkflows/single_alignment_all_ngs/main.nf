

include { SURVIVOR } from '../../modules/survivor/'
include { PhenoSV } from '../../modules/phenosv/'
include { ANNOVAR } from '../../modules/annovar/'
include { Phen2gene } from '../../modules/phen2gene/'
include { RankVar } from '../../modules/rankvar/'
include { Manta } from '../../modules/manta/'
include { xtea } from '../../modules/xtea/'
include { normalize_shortread_alignment } from '../../modules/normalize_shortread_alignment/'
include { CNVnator } from '../../modules/cnvnator/'
include { truvari_shortread_sv_merge } from '../../modules/truvari_shortread_sv_merge/'
include { deepvariant } from '../../modules/deepvariant/'
include { haplotypecaller } from '../../modules/haplotypecaller/'
include { ANNOVAR_SV } from '../../modules/annovar_sv/'
include { common_sv_filter } from '../../modules/common_sv_filter/'
include { ExpansionHunter } from '../../modules/expansion_hunter/'
include { Rankscore_analysis } from '../../modules/rankscore_analysis/'
include { phenotagger } from '../../modules/phenotagger/'
include { phenogpt2 } from '../../modules/phenogpt2/'
include { eh_filter } from '../../modules/eh_filter/'
include { phen2gene_filter } from '../../modules/reduce_region_phen2gene/'
include { ngs_prio } from '../../modules/ngs_prio/'
include { variant_html_report; variant_html_report_with_mito } from '../../modules/variant_html_report/'


// Single sample: short-read full path (SNP + SV + STR) with DeepVariant and Manta.
workflow SINGLE_ALIGNMENT_ALL_NGS {
	take:
	bam
	out_prefix
	ref_fa
	gatk_ref_fa
	eh_variant_catalog
	note
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
	hpo=note
	if ( is_note == "yes" ) {
		if ( params.phenotype_extractor.toString().trim().toLowerCase() == "phenogpt2" ) {
			phenogpt2(note,out_prefix)
			hpo=phenogpt2.out
		}
		else {
			phenotagger(note,out_prefix)
			hpo=phenotagger.out
		}
	}
	Phen2gene(hpo,out_prefix)
	if ( target == "yes" ) {
		phen2_gene_bed=phen2gene_filter(Phen2gene.out,ref_fa,out_prefix,phen2gene_top_n)
		if ( caller_mode == "haplotypecaller" ) {
			haplotypecaller(bam,out_prefix,gatk_ref_fa,phen2_gene_bed)
		}
		else {
			deepvariant(bam,out_prefix,ref_fa,phen2_gene_bed)
		}
	}
	else {
		if ( caller_mode == "haplotypecaller" ) {
			haplotypecaller(bam,out_prefix,gatk_ref_fa,target)
		}
		else {
			deepvariant(bam,out_prefix,ref_fa,target)
		}
	}
	snp_vcf = (caller_mode == "haplotypecaller") ? haplotypecaller.out : deepvariant.out
	annovar_bed = (target == "yes") ? phen2_gene_bed : target
	ANNOVAR(snp_vcf,out_prefix,annovar_bed)
	RankVar(ANNOVAR.out.txt_output,Phen2gene.out,hpo,out_prefix,gnomad,gq,ad,rankvar_filter)
	rankscore_result=Rankscore_analysis(ANNOVAR.out.txt_output,Phen2gene.out,out_prefix,gnomad,rankscore_filter,rankscore_softwares,gq,ad,phen2gene_top_n)
	Manta(bam,out_prefix,ref_fa)
	xtea_mode = params.xtea ? params.xtea.toString().trim().toLowerCase() : "no"
	xtea_vcf = null
	if ( xtea_mode == "yes" ) {
		xtea_input = out_prefix.combine(bam).map { prefix, bam_tuple ->
			tuple([id: prefix], bam_tuple[0], bam_tuple[1])
		}
		xtea(xtea_input, ref_fa)
		xtea_vcf = xtea.out.vcf.map { meta, vcf -> vcf }
	}

	cnvnator_mode = params.cnvnator ? params.cnvnator.toString().trim().toLowerCase() : "yes"
	if ( cnvnator_mode != "no" && xtea_mode == "yes" ) {
		normalize_shortread_alignment(bam,out_prefix,ref_fa)
		CNVnator(normalize_shortread_alignment.out,out_prefix,ref_fa,params.cnvnator_bin_size)
		sv_merge_inputs = Manta.out.combine(CNVnator.out.vcf).combine(xtea_vcf).map { combined_vcfs ->
			combined_vcfs.flatten()
		}
		truvari_shortread_sv_merge(sv_merge_inputs,ref_fa,out_prefix)
		sv_vcf=truvari_shortread_sv_merge.out.merged_vcf
	}
	else if ( cnvnator_mode != "no" ) {
		normalize_shortread_alignment(bam,out_prefix,ref_fa)
		CNVnator(normalize_shortread_alignment.out,out_prefix,ref_fa,params.cnvnator_bin_size)
		sv_merge_inputs = Manta.out.combine(CNVnator.out.vcf).map { combined_vcfs ->
			combined_vcfs.flatten()
		}
		truvari_shortread_sv_merge(sv_merge_inputs,ref_fa,out_prefix)
		sv_vcf=truvari_shortread_sv_merge.out.merged_vcf
	}
	else if ( xtea_mode == "yes" ) {
		sv_merge_inputs = Manta.out.combine(xtea_vcf).map { combined_vcfs ->
			combined_vcfs.flatten()
		}
		truvari_shortread_sv_merge(sv_merge_inputs,ref_fa,out_prefix)
		sv_vcf=truvari_shortread_sv_merge.out.merged_vcf
	}
	else {
		sv_vcf=Manta.out
	}

	sv_annovar_bed = (target == "yes") ? phen2_gene_bed : target
	ANNOVAR_SV(sv_vcf,out_prefix,Phen2gene.out,sv_annovar_bed,"called")
	annovar_sv_for_downstream = ANNOVAR_SV.out
	if ( params.common_sv_filter.toString().trim().toLowerCase() == "yes" ) {
		common_sv_filter(ANNOVAR_SV.out,out_prefix)
		annovar_sv_for_downstream = common_sv_filter.out.filtered_vcf
	}
	SURVIVOR(annovar_sv_for_downstream,out_prefix)
	PhenoSV(SURVIVOR.out.phenosv_inputs,out_prefix,hpo)
	ExpansionHunter(bam,out_prefix,ref_fa,eh_variant_catalog)
	eh_filter(out_prefix,ExpansionHunter.out)
	ngs_prio(out_prefix,RankVar.out,Rankscore_analysis.out.rankscore,Rankscore_analysis.out.clinvar,PhenoSV.out,annovar_sv_for_downstream,ANNOVAR.out.vcf_output,hpo,inheritance_mode,include_clinvar_report,allow_unphased_comphet)
	if ( mito_mode == "yes" ) {
		variant_html_report_with_mito(out_prefix, ngs_prio.out[0], ngs_prio.out[1], eh_filter.out, mito_tsv)
	}
	else {
		variant_html_report(out_prefix, ngs_prio.out[0], ngs_prio.out[1], eh_filter.out)
	}

}
