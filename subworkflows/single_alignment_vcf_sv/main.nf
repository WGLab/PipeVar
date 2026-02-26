

include { SURVIVOR } from '../../modules/survivor/'
include { PhenoSV } from '../../modules/phenosv/'
include { Phen2gene } from '../../modules/phen2gene/'
include { ANNOVAR_SV } from '../../modules/annovar_sv/'
include { phenotagger } from '../../modules/phenotagger/'
include { sv_prio } from '../../modules/sv_prio/'


// Single sample: VCF SV re-annotation and SV prioritization path.
workflow SINGLE_ALIGNMENT_VCF_SV {
	take:
	vcf
	out_prefix
	ref_fa
	note
	is_note
	inheritance_mode

	main:
	
		hpo=note
		if ( is_note == "yes" ) {
			phenotagger(note,out_prefix)
			hpo=phenotagger.out
		}
		Phen2gene(hpo,out_prefix)
		ANNOVAR_SV(vcf,out_prefix,Phen2gene.out)
		SURVIVOR(ANNOVAR_SV.out,out_prefix)
		PhenoSV(SURVIVOR.out,out_prefix,hpo)
		sv_prio(out_prefix,PhenoSV.out,ANNOVAR_SV.out,hpo,inheritance_mode)
}

