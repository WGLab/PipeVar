

include { sniffles } from '../../modules/sniffles/'
include { clair3 } from '../../modules/clair3/'
include { nanocaller } from '../../modules/nanocaller/'
include { CNVpytor } from '../../modules/cnvpytor/'
include { merge_longread_sv_callers } from '../../modules/merge_longread_sv_callers/'
include { SURVIVOR } from '../../modules/survivor/'
include { PhenoSV } from '../../modules/phenosv/'
include { ANNOVAR } from '../../modules/annovar/'
include { Phen2gene } from '../../modules/phen2gene/'
include { RankVar } from '../../modules/rankvar/'
include { NanoRepeat } from '../../modules/nanorepeat/'
include { ANNOVAR_SV } from '../../modules/annovar_sv/'
include { Rankscore_analysis } from '../../modules/rankscore_analysis/'
include { phenotagger } from '../../modules/phenotagger/'
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
		phenotagger(note,out_prefix)
		hpo=phenotagger.out
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
	RankVar(ANNOVAR.out.txt_output,Phen2gene.out,hpo,out_prefix,gnomad,gq,ad,rankvar_filter)
	rankscore_result=Rankscore_analysis(ANNOVAR.out.txt_output,Phen2gene.out,out_prefix,gnomad,rankscore_filter,rankscore_softwares,gq,phen2gene_top_n)
	sniffles(bam,out_prefix,ref_fa)
	cnvpytor_mode = params.cnvpytor ? params.cnvpytor.toString().trim().toLowerCase() : "no"
	cnvpytor_baf_mode = params.cnvpytor_baf ? params.cnvpytor_baf.toString().trim().toLowerCase() : "yes"
	if ( cnvpytor_mode == "yes" ) {
		CNVpytor(bam,out_prefix,snp_vcf,Channel.value(cnvpytor_baf_mode),Channel.value(params.cnvpytor_bin_sizes),Channel.value(params.cnvpytor_primary_bin),Channel.value(params.cnvpytor_min_size))
		sv_merge_inputs = sniffles.out.combine(CNVpytor.out.vcf).map { combined_vcfs ->
			combined_vcfs.flatten()
		}
		merge_longread_sv_callers(sv_merge_inputs,out_prefix)
		sv_vcf = merge_longread_sv_callers.out
	}
	else {
		sv_vcf = sniffles.out
	}
	annovar_sv_bed = (target == "yes") ? phen2_gene_bed : target
	ANNOVAR_SV(sv_vcf,out_prefix,Phen2gene.out,annovar_sv_bed)
	SURVIVOR(ANNOVAR_SV.out,out_prefix)
	PhenoSV(SURVIVOR.out,out_prefix,hpo)
	NanoRepeat(bam,out_prefix,ref_fa)
	longphase(bam,ANNOVAR.out.vcf_output,ANNOVAR_SV.out,PhenoSV.out,rankscore_result.rankscore,rankscore_result.clinvar,RankVar.out,hpo,out_prefix,ref_fa,inheritance_mode,include_clinvar_report,allow_unphased_comphet)
	if ( mito_mode == "yes" ) {
		variant_html_report_with_mito(out_prefix, longphase.out[0], longphase.out[1], NanoRepeat.out, mito_tsv)
	}
	else {
		variant_html_report(out_prefix, longphase.out[0], longphase.out[1], NanoRepeat.out)
	}

}
