

include { nanocaller } from '../../modules/nanocaller/'
include { ANNOVAR } from '../../modules/annovar/'
include { Phen2gene } from '../../modules/phen2gene/'
include { RankVar } from '../../modules/rankvar/'
include { Rankscore_analysis } from '../../modules/rankscore_analysis/'
include { phenotagger } from '../../modules/phenotagger/'
include { phen2gene_filter } from '../../modules/reduce_region_phen2gene/'
include { snp_prio } from '../../modules/snp_prio/'



// Single sample: long-read SNP-only light path using NanoCaller.
workflow SINGLE_ALIGNMENT_LONG_SNP_LIGHT {
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

	rankvar_filter
	is_note
	target
	inheritance_mode

	main:
	
	hpo=note
	if ( is_note == "yes" ) {
		phenotagger(note,out_prefix)
		hpo=phenotagger.out
	}
	Phen2gene(hpo,out_prefix)
        if ( target == "yes" ) {
                phen2_gene_bed=phen2gene_filter(Phen2gene.out,ref_fa,out_prefix,phen2gene_top_n)
                nanocaller(bam,out_prefix,ref_fa,phen2_gene_bed)
        }
        else {
                nanocaller(bam,out_prefix,ref_fa,target)
        }
	annovar_bed = (target == "yes") ? phen2_gene_bed : target
	ANNOVAR(nanocaller.out,out_prefix,annovar_bed)
	RankVar(ANNOVAR.out.txt_output,Phen2gene.out,hpo,out_prefix,gnomad,gq,ad,rankvar_filter)
	Rankscore_analysis(ANNOVAR.out.txt_output,Phen2gene.out,out_prefix,gnomad,rankscore_filter,gq,phen2gene_top_n)
	snp_prio(out_prefix,Rankscore_analysis.out.rankscore,Rankscore_analysis.out.clinvar,RankVar.out,ANNOVAR.out.vcf_output,hpo,inheritance_mode)
}
