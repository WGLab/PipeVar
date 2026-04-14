include { mito_prep_mutect2 } from '../../modules/mito_prep_mutect2/'
include { mito_mutect2 } from '../../modules/mito_mutect2/'
include { mito_annotation } from '../../modules/mito_annotation/'
include { mito_prio } from '../../modules/mito_prio/'

workflow SINGLE_ALIGNMENT_NGS_MITO {
	take:
	bam
	out_prefix
	ref_fa
	mito_contig

	main:
	mutect2_ref_fa = ref_fa.map { fa_file, fai_file, dict_file, bwa_amb, bwa_ann, bwa_bwt, bwa_pac, bwa_sa ->
		tuple(fa_file, fai_file, dict_file)
	}
	mito_prep_mutect2(bam, out_prefix, ref_fa, mito_contig)
	mito_mutect2(mito_prep_mutect2.out, out_prefix, mutect2_ref_fa, mito_contig)
	mito_annotation(mito_mutect2.out[0], out_prefix)
	mito_prio(mito_annotation.out[0], out_prefix)

	emit:
	mito_vcf = mito_mutect2.out[0]
	annotated_tsv = mito_annotation.out[0]
	prioritized_tsv = mito_prio.out
}
