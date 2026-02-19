

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
include { sv_prio } from '../../modules/sv_prio/'


workflow SINGLE_ALIGNMENT_NGS_SV {
	take:
	bam
	out_prefix
	ref_fa
	note
	is_note

	main:

	hpo=note
	if ( is_note == "yes" ) {
		phenotagger(note,out_prefix)
		hpo=phenotagger.out
	}
	Manta(bam,out_prefix,ref_fa)
	ANNOVAR_SV(out_prefix,Manta.out)
	SURVIVOR(ANNOVAR_SV.out,out_prefix)
	PhenoSV(SURVIVOR.out,out_prefix,hpo)
	ExpansionHunter(bam,out_prefix,ref_fa)
	eh_filter(out_prefix,ExpansionHunter.out)
	sv_prio(out_prefix,ANNOVAR_SV.out,PhenoSV.out)
}



