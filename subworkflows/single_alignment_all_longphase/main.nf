

include { cuteSV } from '../../modules/cutesv/'
include { sniffles } from '../../modules/sniffles/'
include { nanocaller } from '../../modules/nanocaller/'
include { clair3 } from '../../modules/clair3/'
include { SURVIVOR } from '../../modules/survivor/'
include { PhenoSV } from '../../modules/phenosv/'
include { ANNOVAR } from '../../modules/annovar/'
include { Phen2gene } from '../../modules/phen2gene/'
include { RankVar } from '../../modules/rankvar/'
include { NanoRepeat } from '../../modules/nanorepeat/'
include { truvari } from '../../modules/truvari/'
include { Manta } from '../../modules/manta/'
include { haplotypecaller } from '../../modules/haplotypecaller/'
include { deepvariant } from '../../modules/deepvariant/'
include { ANNOVAR_SV } from '../../modules/annovar_sv/'
include { ExpansionHunter } from '../../modules/expansion_hunter/'
include { Rankscore_analysis } from '../../modules/rankscore_analysis/'
include { phenotagger } from '../../modules/phenotagger/'
include { longphase } from '../../modules/longphase/'
include { eh_filter } from '../../modules/eh_filter/'
include { phen2gene_filter } from '../../modules/reduce_region_phen2gene/'

workflow SINGLE_ALIGNMENT_ALL_LONGPHASE {
	take:
	bam
	out_prefix
	ref_fa
	note
	rankscore_filter
	gnomad
	gq
	ad
	is_note
	target

	main:
	
	hpo=note
	if ( is_note == "yes" ) {
		hpo=phenotagger(note,out_prefix)
	}
	Phen2gene(hpo,out_prefix)
	if ( target == "yes" ) {
		phen2_gene_bed=phen2gene_filter(Phen2gene.out,ref_fa,out_prefix)
		clair3(bam,out_prefix,ref_fa,phen2_gene_bed)
	}
	else {
		clair3(bam,out_prefix,ref_fa,target)
	}
	ANNOVAR(clair3.out,out_prefix)
	RankVar(ANNOVAR.out.txt_output,Phen2gene.out,hpo,out_prefix,gnomad,gq,ad)
	rankscore_result=Rankscore_analysis(ANNOVAR.out.txt_output,Phen2gene.out,out_prefix,gnomad,rankscore_filter)
	cuteSV(bam,out_prefix,ref_fa)
	sniffles(bam,out_prefix,ref_fa)
	truvari(cuteSV.out,sniffles.out,ref_fa,out_prefix)
	SURVIVOR(truvari.out,out_prefix)
	PhenoSV(SURVIVOR.out,out_prefix,hpo)
	NanoRepeat(bam,out_prefix,ref_fa)
	longphase(bam,ANNOVAR.out.vcf_output,SURVIVOR.out,PhenoSV.out,rankscore_result.rankscore,rankscore_result.clinvar,RankVar.out,out_prefix,ref_fa)

}



