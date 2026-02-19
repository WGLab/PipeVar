

include { ANNOVAR } from '../../modules/annovar/'
include { Phen2gene } from '../../modules/phen2gene/'
include { RankVar } from '../../modules/rankvar/'
include { haplotypecaller } from '../../modules/haplotypecaller/'
include { ExpansionHunter } from '../../modules/expansion_hunter/'
include { Rankscore_analysis } from '../../modules/rankscore_analysis/'
include { phenotagger } from '../../modules/phenotagger/'
include { eh_filter } from '../../modules/eh_filter/'
include { phen2gene_filter } from '../../modules/reduce_region_phen2gene/'
include { snp_prio } from '../../modules/snp_prio/'


// Single sample: short-read SNP-only light path with HaplotypeCaller.
workflow SINGLE_ALIGNMENT_NGS_SNP_LIGHT {
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

	main:

	hpo=note
	if ( is_note == "yes" ) {
		phenotagger(note,out_prefix)
		hpo=phenotagger.out
	}
	Phen2gene(hpo,out_prefix)
	if ( target == "yes" ) {
		phen2_gene_bed=phen2gene_filter(Phen2gene.out,ref_fa,out_prefix,phen2gene_top_n)
		haplotypecaller(bam,out_prefix,ref_fa,phen2_gene_bed)
	}
	else {
		haplotypecaller(bam,out_prefix,ref_fa,target)
	}
	ANNOVAR(haplotypecaller.out,out_prefix)
	RankVar(ANNOVAR.out.txt_output,Phen2gene.out,hpo,out_prefix,gnomad,gq,ad)
	Rankscore_analysis(ANNOVAR.out.txt_output,Phen2gene.out,out_prefix,gnomad,rankscore_filter,phen2gene_top_n)
	ExpansionHunter(bam,out_prefix,ref_fa)
	eh_filter(out_prefix,ExpansionHunter.out)
	snp_prio(out_prefix,Rankscore_analysis.out.rankscore,Rankscore_analysis.out.clinvar,RankVar.out,ANNOVAR.out.vcf_output)
}


