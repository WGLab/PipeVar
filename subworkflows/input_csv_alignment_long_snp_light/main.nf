
include { multi_annovar } from '../../modules/multi_annovar/'
include { multi_phen2gene } from '../../modules/multi_phen2gene/'
include { multi_rankscore } from '../../modules/multi_rankscore/'
include { multi_annovar_sv } from '../../modules/multi_annovar_sv/'
include { multi_survivor } from '../../modules/multi_survivor/'
include { multi_phenosv } from '../../modules/multi_phenosv/'
include { multi_rankvar } from '../../modules/multi_rankvar/'
include { multi_sniffles } from '../../modules/multi_sniffles/'
include { multi_manta } from '../../modules/multi_manta/'
include { multi_clair3 } from '../../modules/multi_clair3/'
include { multi_nanorepeat } from '../../modules/multi_nanorepeat/'
include { multi_expansionhunter } from '../../modules/multi_expansionhunter/'
include { multi_eh_filter } from '../../modules/multi_eh_filter/'
include { multi_deepvariant } from '../../modules/multi_deepvariant/'
include { multi_nanocaller } from '../../modules/multi_nanocaller/'
include { multi_haplotypecaller } from '../../modules/multi_haplotypecaller/'
include { multi_prep_gatk } from '../../modules/multi_prep_gatk/'
include { multi_longphase } from '../../modules/multi_longphase/'
include { multi_phenotagger } from '../../modules/multi_phenotagger/'
include { multi_phen2gene_filter } from '../../modules/multi_reduce_region_phen2gene/'
include { multi_snp_prio } from '../../modules/multi_snp_prio/'



workflow INPUT_CSV_ALIGNMENT_LONG_SNP_LIGHT {
	take:
	input_bam
	ref_fa
	rankscore_filter
	gnomad
	gq
	ad
	is_note
	target

	main:

        input_bam_no_bam =  input_bam.map { out_prefix, bam_file, bai_file, note_file -> return tuple ( out_prefix,note_file ) }
        input_bam_with_bam= input_bam.map { out_prefix, bam_file, bai_file, note_file -> return tuple (out_prefix, bam_file, bai_file) }
        if ( is_note == "yes" ) {
                input_bam_no_bam=multi_phenotagger(input_bam_no_bam)
        }
	phen2gene_result=multi_phen2gene(input_bam_no_bam)
        if ( target == "yes" ) {
                phen2_gene_bed=multi_phen2gene_filter(phen2gene_result,ref_fa)
                nanocaller_result=multi_nanocaller(input_bam_with_bam,ref_fa,phen2_gene_bed)
        }
        else {
                nanocaller_result=multi_nanocaller(input_bam_with_bam,ref_fa,target)
        }
	annovar_result=multi_annovar(nanocaller_result)
	annovar_result_txt=annovar_result.map { item -> tuple(item[0], item[1]) }
	join_annovar_phen2gene=annovar_result_txt.join(phen2gene_result)
	join_annovar_hpo=join_annovar_phen2gene.join(input_bam_no_bam)
	rankscore_result=multi_rankscore(join_annovar_phen2gene,gnomad,rankscore_filter)
	rankvar_result=multi_rankvar(join_annovar_hpo,gnomad,gq,ad)
	rankscore_rankvar_join=rankvar_result.join(rankscore_result)
        annovar_result_vcf=annovar_result.map { item -> tuple(item[0], item[2]) }
        snp_prio_input=rankscore_rankvar_join.join(annovar_result_vcf)
        multi_snp_prio(snp_prio_input)

}	



