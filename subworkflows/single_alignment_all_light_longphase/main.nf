

include { sniffles } from '../../modules/sniffles/'
include { nanocaller } from '../../modules/nanocaller/'
include { SURVIVOR } from '../../modules/survivor/'
include { PhenoSV } from '../../modules/phenosv/'
include { ANNOVAR } from '../../modules/annovar/'
include { Phen2gene } from '../../modules/phen2gene/'
include { RankVar } from '../../modules/rankvar/'
include { NanoRepeat } from '../../modules/nanorepeat/'
include { ANNOVAR_SV } from '../../modules/annovar_sv/'
include { Rankscore_analysis } from '../../modules/rankscore_analysis/'
include { phenotagger } from '../../modules/phenotagger/'
include { phenogpt2 } from '../../modules/phenogpt2/'
include { longphase } from '../../modules/longphase/'
include { phen2gene_filter } from '../../modules/reduce_region_phen2gene/'
include { variant_html_report } from '../../modules/variant_html_report/'


// Single sample: light long-read full path (SNP + SV + STR) with NanoCaller and longphase.
workflow SINGLE_ALIGNMENT_ALL_LIGHT_LONGPHASE {
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
                nanocaller(bam,out_prefix,ref_fa,phen2_gene_bed)
        }
        else {
                nanocaller(bam,out_prefix,ref_fa,target)
        }
	annovar_bed = (target == "yes") ? phen2_gene_bed : target
	ANNOVAR(nanocaller.out,out_prefix,annovar_bed)
	RankVar(ANNOVAR.out.txt_output,Phen2gene.out,hpo,out_prefix,gnomad,gq,ad,rankvar_filter)
	rankscore_result=Rankscore_analysis(ANNOVAR.out.txt_output,Phen2gene.out,out_prefix,gnomad,rankscore_filter,rankscore_softwares,gq,phen2gene_top_n)
	sniffles(bam,out_prefix,ref_fa)
	annovar_sv_bed = (target == "yes") ? phen2_gene_bed : target
	ANNOVAR_SV(sniffles.out,out_prefix,Phen2gene.out,annovar_sv_bed,"called")
	SURVIVOR(ANNOVAR_SV.out,out_prefix)
	PhenoSV(SURVIVOR.out,out_prefix,hpo)
	NanoRepeat(bam,out_prefix,ref_fa)
	longphase(bam,ANNOVAR.out.vcf_output,sniffles.out,ANNOVAR_SV.out,PhenoSV.out,rankscore_result.rankscore,rankscore_result.clinvar,RankVar.out,hpo,out_prefix,ref_fa,inheritance_mode,include_clinvar_report,allow_unphased_comphet)
	variant_html_report(out_prefix, longphase.out[0], longphase.out[1], NanoRepeat.out)

}
