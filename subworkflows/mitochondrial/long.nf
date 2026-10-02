include { mito_clair3 } from '../../modules/mito_clair3'
include { mito_clair3_postprocess } from '../../modules/mito_clair3_postprocess'
include { mito_annotation } from '../../modules/mito_annotation'
include { mito_prio } from '../../modules/mito_prio'

// Sample-keyed long-read mitochondrial calling, annotation, and prioritization.
workflow ALIGNMENT_LONG_MITO {
    take:
    alignments
    reference_bundle
    mito_contig

    main:
    raw_vcf = mito_clair3(alignments, reference_bundle, mito_contig)
    mito_vcf = mito_clair3_postprocess(raw_vcf)
    annotated = mito_annotation(mito_vcf)
    prioritized = mito_prio(annotated)

    emit:
    mito_vcf = mito_vcf
    annotated_tsv = annotated
    prioritized_tsv = prioritized
}
