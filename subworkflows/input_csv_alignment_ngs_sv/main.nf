
include { multi_phen2gene } from '../../modules/multi_phen2gene/'
include { multi_annovar_sv } from '../../modules/multi_annovar_sv/'
include { multi_survivor } from '../../modules/multi_survivor/'
include { multi_phenosv } from '../../modules/multi_phenosv/'
include { multi_manta } from '../../modules/multi_manta/'
include { multi_expansionhunter } from '../../modules/multi_expansionhunter/'
include { multi_eh_filter } from '../../modules/multi_eh_filter/'
include { multi_phenotagger } from '../../modules/multi_phenotagger/'
include { multi_sv_prio } from '../../modules/multi_sv_prio/'



// CSV batch: short-read SV/STR prioritization path with Manta.
workflow INPUT_CSV_ALIGNMENT_NGS_SV {
	take:
	input_bam
	ref_fa
	is_note

	main:

        input_bam_no_bam =  input_bam.map { out_prefix, bam_file, bai_file, note_file -> return tuple ( out_prefix,note_file ) }
        input_bam_with_bam= input_bam.map { out_prefix, bam_file, bai_file, note_file -> return tuple (out_prefix, bam_file, bai_file) }
        if ( is_note == "yes" ) {
                input_bam_no_bam=multi_phenotagger(input_bam_no_bam)
        }
        phen2gene_result=multi_phen2gene(input_bam_no_bam)
        multi_eh_result=multi_expansionhunter(input_bam_with_bam,ref_fa)
        multi_eh_filter(multi_eh_result)
        manta_result=multi_manta(input_bam_with_bam,ref_fa)
        manta_result_annovar=manta_result.join(input_bam_no_bam)
	annovar_sv_result=multi_annovar_sv(manta_result_annovar)
	survivor_result=multi_survivor(annovar_sv_result)
        phenosv_input=survivor_result.join(input_bam_no_bam)
        phenosv_result=multi_phenosv(phenosv_input)
	sv_prio_input=phenosv_result.join(annovar_sv_result)
        multi_sv_prio(sv_prio_input)

}	



