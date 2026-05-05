include { multi_mito_clair3 } from '../../modules/multi_mito_clair3/'
include { multi_mito_annotation } from '../../modules/multi_mito_annotation/'
include { multi_mito_prio } from '../../modules/multi_mito_prio/'

workflow INPUT_CSV_ALIGNMENT_LONG_MITO {
	take:
	input_bam
	ref_fa
	mito_contig

	main:
	input_bam_with_bam = input_bam.map { out_prefix, bam_file, bai_file, note_file -> tuple(out_prefix, bam_file, bai_file) }
	mt_vcf = multi_mito_clair3(input_bam_with_bam, ref_fa, mito_contig)
	annotated = multi_mito_annotation(mt_vcf)
	multi_mito_prio(annotated)

	emit:
	mito_vcf = mt_vcf
	annotated_tsv = annotated
	prioritized_tsv = multi_mito_prio.out
}
