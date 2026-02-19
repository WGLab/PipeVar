

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


workflow SINGLE_ALIGNMENT_ALL_NGS {
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
		phenotagger(note,out_prefix)
		hpo=phenotagger.out
	}
	Phen2gene(hpo,out_prefix)
	if ( target == "yes" ) {
		phen2_gene_bed=phen2gene_filter(Phen2gene.out,ref_fa,out_prefix)
		deepvariant(bam,out_prefix,ref_fa,phen2_gene_bed)
	}
	else {
		deepvariant(bam,out_prefix,ref_fa,target)
	}
	ANNOVAR(deepvariant.out,out_prefix)
	RankVar(ANNOVAR.out.txt_output,Phen2gene.out,hpo,out_prefix,gnomad,gq,ad)
	rankscore_result=Rankscore_analysis(ANNOVAR.out.txt_output,Phen2gene.out,out_prefix,gnomad,rankscore_filter)
	Manta(bam,out_prefix,ref_fa)
	SURVIVOR(Manta.out,out_prefix)
	PhenoSV(SURVIVOR.out,out_prefix,hpo)
	ExpansionHunter(bam,out_prefix,ref_fa)
	eh_filter(out_prefix,ExpansionHunter.out)

}



