include { multi_mito_prep_mutect2 } from '../../modules/multi_mito_prep_mutect2/'
include { multi_mito_mutect2 } from '../../modules/multi_mito_mutect2/'
include { multi_mito_annotation } from '../../modules/multi_mito_annotation/'
include { multi_mito_prio } from '../../modules/multi_mito_prio/'

workflow INPUT_CSV_ALIGNMENT_NGS_MITO {
	take:
	input_bam
	ref_fa
	mito_contig

	main:
	input_bam_with_bam = input_bam.map { out_prefix, bam_file, bai_file, note_file -> tuple(out_prefix, bam_file, bai_file) }
	mutect2_ref_fa = ref_fa.map { fa_file, fai_file, dict_file, bwa_amb, bwa_ann, bwa_bwt, bwa_pac, bwa_sa ->
		tuple(fa_file, fai_file, dict_file)
	}
	prepped = multi_mito_prep_mutect2(input_bam_with_bam, ref_fa, mito_contig)
	mt_vcf = multi_mito_mutect2(prepped, mutect2_ref_fa, mito_contig)
	annotated = multi_mito_annotation(mt_vcf)
	multi_mito_prio(annotated)

	emit:
	mito_vcf = mt_vcf
	annotated_tsv = annotated
	prioritized_tsv = multi_mito_prio.out
}
