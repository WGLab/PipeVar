

include { ANNOVAR } from '../../modules/annovar/'
include { Phen2gene } from '../../modules/phen2gene/'
include { RankVar } from '../../modules/rankvar/'
include { Rankscore_analysis } from '../../modules/rankscore_analysis/'
include { phenotagger } from '../../modules/phenotagger/'
include { snp_prio } from '../../modules/snp_prio/'

// Single sample: VCF SNP re-annotation and SNP prioritization path.
workflow SINGLE_ALIGNMENT_VCF_SNP {
	take:
	vcf
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
	inheritance_mode

	main:
	
	hpo=note
	if ( is_note == "yes" ) {
		phenotagger(note,out_prefix)
		hpo=phenotagger.out
	}
	ANNOVAR(vcf,out_prefix)
	Phen2gene(hpo,out_prefix)
	RankVar(ANNOVAR.out.txt_output,Phen2gene.out,hpo,out_prefix,gnomad,gq,ad,rankvar_filter)
	Rankscore_analysis(ANNOVAR.out.txt_output,Phen2gene.out,out_prefix,gnomad,rankscore_filter,gq,phen2gene_top_n)
	snp_prio(out_prefix,Rankscore_analysis.out.rankscore,Rankscore_analysis.out.clinvar,RankVar.out,ANNOVAR.out.vcf_output,hpo,inheritance_mode)
}


