

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
include { snp_prio } from '../../modules/snp_prio/'

workflow SINGLE_ALIGNMENT_VCF_SNP {
	take:
	vcf
	out_prefix
	ref_fa
	note
	rankscore_filter
	gnomad
	gq
	ad
	is_note

	main:
	
	hpo=note
	if ( is_note == "yes" ) {
		phenotagger(note,out_prefix)
		hpo=phenotagger.out
	}
	ANNOVAR(vcf,out_prefix)
	Phen2gene(hpo,out_prefix)
	RankVar(ANNOVAR.out.txt_output,Phen2gene.out,hpo,out_prefix,gnomad,gq,ad)
	Rankscore_analysis(ANNOVAR.out.txt_output,Phen2gene.out,out_prefix,gnomad,rankscore_filter)
	snp_prio(out_prefix,Rankscore_analysis.out.rankscore,Rankscore_analysis.out.clinvar,RankVar.out,ANNOVAR.out.vcf_output)
}



