

include { sniffles } from '../../modules/sniffles/'
include { clair3 } from '../../modules/clair3/'
include { nanocaller } from '../../modules/nanocaller/'
include { SURVIVOR } from '../../modules/survivor/'
include { PhenoSV } from '../../modules/phenosv/'
include { ANNOVAR } from '../../modules/annovar/'
include { Phen2gene } from '../../modules/phen2gene/'
include { RankVar } from '../../modules/rankvar/'
include { RankVar as RankVarNanoCaller } from '../../modules/rankvar/'
include { NanoRepeat } from '../../modules/nanorepeat/'
include { ANNOVAR_SV } from '../../modules/annovar_sv/'
include { common_sv_filter } from '../../modules/common_sv_filter/'
include { Rankscore_analysis } from '../../modules/rankscore_analysis/'
include { phenotagger } from '../../modules/phenotagger/'
include { phenogpt2 } from '../../modules/phenogpt2/'
include { longphase } from '../../modules/longphase/'
include { phen2gene_filter } from '../../modules/reduce_region_phen2gene/'
include { variant_html_report; variant_html_report_with_mito } from '../../modules/variant_html_report/'

// Single sample: long-read full path (SNP + SV + STR) with Clair3/Sniffles and longphase.
workflow SINGLE_ALIGNMENT_ALL_LONGPHASE {
	take:
	bam
	out_prefix
	ref_fa
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
		if ( caller_mode == "nanocaller" ) {
			nanocaller(bam,out_prefix,ref_fa,phen2_gene_bed)
		}
		else {
			clair3(bam,out_prefix,ref_fa,phen2_gene_bed)
		}
	}
	else {
		if ( caller_mode == "nanocaller" ) {
			nanocaller(bam,out_prefix,ref_fa,target)
		}
		else {
			clair3(bam,out_prefix,ref_fa,target)
		}
	}
	snp_vcf = (caller_mode == "nanocaller") ? nanocaller.out : clair3.out
	annovar_bed = (target == "yes") ? phen2_gene_bed : target
	ANNOVAR(snp_vcf,out_prefix,annovar_bed)
	if ( caller_mode == "nanocaller" ) {
		RankVarNanoCaller(ANNOVAR.out.txt_output,Phen2gene.out,hpo,out_prefix,gnomad,gq,ad,rankvar_filter)
		rankvar_result = RankVarNanoCaller.out
	}
	else {
		RankVar(ANNOVAR.out.txt_output,Phen2gene.out,hpo,out_prefix,gnomad,gq,ad,rankvar_filter)
		rankvar_result = RankVar.out
	}
	rankscore_result=Rankscore_analysis(ANNOVAR.out.txt_output,Phen2gene.out,out_prefix,gnomad,rankscore_filter,rankscore_softwares,gq,ad,phen2gene_top_n)
	sniffles(bam,out_prefix,ref_fa)
	ANNOVAR_SV(sniffles.out,out_prefix,"called")
	annovar_sv_for_downstream = ANNOVAR_SV.out
	if ( params.common_sv_filter.toString().trim().toLowerCase() == "yes" ) {
		common_sv_filter(ANNOVAR_SV.out,out_prefix)
		annovar_sv_for_downstream = common_sv_filter.out.filtered_vcf
	}
	SURVIVOR(annovar_sv_for_downstream,out_prefix)
	PhenoSV(SURVIVOR.out.phenosv_inputs,out_prefix,hpo)
	NanoRepeat(bam,out_prefix,ref_fa)
	longphase(bam,ANNOVAR.out.vcf_output,sniffles.out,annovar_sv_for_downstream,PhenoSV.out,rankscore_result.rankscore,rankscore_result.clinvar,rankvar_result,hpo,out_prefix,ref_fa,inheritance_mode,include_clinvar_report,allow_unphased_comphet)
	if ( mito_mode == "yes" ) {
		variant_html_report_with_mito(out_prefix, longphase.out[0], longphase.out[1], NanoRepeat.out, mito_tsv)
	}
	else {
		variant_html_report(out_prefix, longphase.out[0], longphase.out[1], NanoRepeat.out)
	}

}
