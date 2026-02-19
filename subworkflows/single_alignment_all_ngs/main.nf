

include { SURVIVOR } from '../../modules/survivor/'
include { PhenoSV } from '../../modules/phenosv/'
include { ANNOVAR } from '../../modules/annovar/'
include { Phen2gene } from '../../modules/phen2gene/'
include { RankVar } from '../../modules/rankvar/'
include { Manta } from '../../modules/manta/'
include { deepvariant } from '../../modules/deepvariant/'
include { haplotypecaller } from '../../modules/haplotypecaller/'
include { ExpansionHunter } from '../../modules/expansion_hunter/'
include { Rankscore_analysis } from '../../modules/rankscore_analysis/'
include { phenotagger } from '../../modules/phenotagger/'
include { eh_filter } from '../../modules/eh_filter/'
include { phen2gene_filter } from '../../modules/reduce_region_phen2gene/'


// Single sample: short-read full path (SNP + SV + STR) with DeepVariant and Manta.
workflow SINGLE_ALIGNMENT_ALL_NGS {
	take:
	bam
	out_prefix
	ref_fa
	note
	rankscore_filter
	phen2gene_top_n
	gnomad
	gq
	ad
	is_note
	target
	caller_mode

	main:

	hpo=note
	if ( is_note == "yes" ) {
		phenotagger(note,out_prefix)
		hpo=phenotagger.out
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
	ANNOVAR(snp_vcf,out_prefix)
	RankVar(ANNOVAR.out.txt_output,Phen2gene.out,hpo,out_prefix,gnomad,gq,ad)
	rankscore_result=Rankscore_analysis(ANNOVAR.out.txt_output,Phen2gene.out,out_prefix,gnomad,rankscore_filter,phen2gene_top_n)
	Manta(bam,out_prefix,ref_fa)
	SURVIVOR(Manta.out,out_prefix)
	PhenoSV(SURVIVOR.out,out_prefix,hpo)
	ExpansionHunter(bam,out_prefix,ref_fa)
	eh_filter(out_prefix,ExpansionHunter.out)

}


