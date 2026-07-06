

include { ANNOVAR } from '../../modules/annovar/'
include { Phen2gene } from '../../modules/phen2gene/'
include { RankVar } from '../../modules/rankvar/'
include { deepvariant } from '../../modules/deepvariant/'
include { haplotypecaller } from '../../modules/haplotypecaller/'
include { ExpansionHunter } from '../../modules/expansion_hunter/'
include { Rankscore_analysis } from '../../modules/rankscore_analysis/'
include { phenotagger } from '../../modules/phenotagger/'
include { phenogpt2 } from '../../modules/phenogpt2/'
include { eh_filter } from '../../modules/eh_filter/'
include { phen2gene_filter } from '../../modules/reduce_region_phen2gene/'
include { snp_prio } from '../../modules/snp_prio/'


// Single sample: short-read SNP-only prioritization path with DeepVariant.
workflow SINGLE_ALIGNMENT_NGS_SNP {
	take:
	bam
	out_prefix
	ref_fa
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
			haplotypecaller(bam,out_prefix,ref_fa,phen2_gene_bed)
		}
		else {
			deepvariant(bam,out_prefix,ref_fa,phen2_gene_bed)
		}
	}
	else {
		if ( caller_mode == "haplotypecaller" ) {
			haplotypecaller(bam,out_prefix,ref_fa,target)
		}
		else {
			deepvariant(bam,out_prefix,ref_fa,target)
		}
	}
	snp_vcf = (caller_mode == "haplotypecaller") ? haplotypecaller.out : deepvariant.out
	annovar_bed = (target == "yes") ? phen2_gene_bed : target
	ANNOVAR(snp_vcf,out_prefix,annovar_bed)
	RankVar(ANNOVAR.out.txt_output,Phen2gene.out,hpo,out_prefix,gnomad,gq,ad,rankvar_filter)
	Rankscore_analysis(ANNOVAR.out.txt_output,Phen2gene.out,out_prefix,gnomad,rankscore_filter,rankscore_softwares,gq,phen2gene_top_n)
	ExpansionHunter(bam,out_prefix,ref_fa,eh_variant_catalog)
	eh_filter(out_prefix,ExpansionHunter.out)
	snp_prio(out_prefix,Rankscore_analysis.out.rankscore,Rankscore_analysis.out.clinvar,RankVar.out,ANNOVAR.out.vcf_output,hpo,inheritance_mode,include_clinvar_report,allow_unphased_comphet)
}
