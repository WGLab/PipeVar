include { mito_prep_mutect2 } from '../../modules/mito_prep_mutect2'
include { mito_mutect2 } from '../../modules/mito_mutect2'
include { mito_annotation } from '../../modules/mito_annotation'
include { mito_prio } from '../../modules/mito_prio'

// Sample-keyed short-read mitochondrial calling, annotation, and prioritization.
workflow ALIGNMENT_NGS_MITO {
    take:
    alignments
    reference_bundle
    mito_contig

    main:
    mutect2_ref_fa = reference_bundle.map { fa_file, fai_file, dict_file, bwa_amb, bwa_ann, bwa_bwt, bwa_pac, bwa_sa ->
        tuple(fa_file, fai_file, dict_file)
    }
    prepped = mito_prep_mutect2(alignments, reference_bundle, mito_contig)
    mito_vcf = mito_mutect2(prepped, mutect2_ref_fa, mito_contig)
    annotated = mito_annotation(mito_vcf)
    prioritized = mito_prio(annotated)

    emit:
    mito_vcf = mito_vcf
    annotated_tsv = annotated
    prioritized_tsv = prioritized
}
