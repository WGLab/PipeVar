include { mito_clair3 } from '../../modules/mito_clair3/'
include { mito_clair3_postprocess } from '../../modules/mito_clair3_postprocess/'
include { mito_annotation } from '../../modules/mito_annotation/'
include { mito_prio } from '../../modules/mito_prio/'

workflow SINGLE_ALIGNMENT_LONG_MITO {
	take:
	bam
	out_prefix
	ref_fa
	mito_contig

	main:
	mito_clair3(bam, out_prefix, ref_fa, mito_contig)
	mito_clair3_postprocess(mito_clair3.out[0], out_prefix)
	mito_annotation(mito_clair3_postprocess.out[0], out_prefix)
	mito_prio(mito_annotation.out[0], out_prefix)

	emit:
	mito_vcf = mito_clair3_postprocess.out[0]
	annotated_tsv = mito_annotation.out[0]
	prioritized_tsv = mito_prio.out
}
