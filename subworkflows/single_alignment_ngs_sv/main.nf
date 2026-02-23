

include { SURVIVOR } from '../../modules/survivor/'
include { PhenoSV } from '../../modules/phenosv/'
include { Phen2gene } from '../../modules/phen2gene/'
include { Manta } from '../../modules/manta/'
include { ANNOVAR_SV } from '../../modules/annovar_sv/'
include { ExpansionHunter } from '../../modules/expansion_hunter/'
include { phenotagger } from '../../modules/phenotagger/'
include { eh_filter } from '../../modules/eh_filter/'
include { sv_prio } from '../../modules/sv_prio/'


// Single sample: short-read SV/STR prioritization path with Manta -> ANNOVAR_SV.
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
	Phen2gene(hpo,out_prefix)
	Manta(bam,out_prefix,ref_fa)
		ANNOVAR_SV(Manta.out,out_prefix,Phen2gene.out)
	SURVIVOR(ANNOVAR_SV.out,out_prefix)
	PhenoSV(SURVIVOR.out,out_prefix,hpo)
	ExpansionHunter(bam,out_prefix,ref_fa)
	eh_filter(out_prefix,ExpansionHunter.out)
		sv_prio(out_prefix,PhenoSV.out,ANNOVAR_SV.out,hpo)
}

