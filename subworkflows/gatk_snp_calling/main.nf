include { GATK_RESOURCE_BUNDLE } from '../../modules/local/gatk_resource_bundle'
include { GATK4_HAPLOTYPECALLER } from '../../modules/nf-core/gatk4/haplotypecaller'
include { GATK4_VARIANTRECALIBRATOR } from '../../modules/nf-core/gatk4/variantrecalibrator'
include { GATK4_APPLYVQSR } from '../../modules/nf-core/gatk4/applyvqsr'

// Call and recalibrate keyed short-read samples with one validated hg38 bundle
// shared across the run. Optional regions are represented by [].
workflow GATK_SMALL_VARIANT_CALLING {
    take:
    // (sample_id, BAM/CRAM, index, regions); [] requests whole-reference calling.
    input_bam_regions
    reference_bundle

    main:
    reference = reference_bundle
    reference_fasta = reference.map { fasta, fai, dict -> tuple([id: 'reference'], fasta) }
    reference_fai = reference.map { fasta, fai, dict -> tuple([id: 'reference'], fai) }
    reference_dict = reference.map { fasta, fai, dict -> tuple([id: 'reference'], dict) }

    haplotypecaller_input = input_bam_regions.map { sample_id, bam, bai, regions ->
        tuple([id: sample_id, single_end: false], bam, bai, regions, [])
    }
    empty_dbsnp = Channel.value(tuple([id: 'reference'], []))
    empty_dbsnp_tbi = Channel.value(tuple([id: 'reference'], []))

    GATK4_HAPLOTYPECALLER(
        haplotypecaller_input,
        reference_fasta,
        reference_fai,
        reference_dict,
        empty_dbsnp,
        empty_dbsnp_tbi
    )

    // One eight-file VQSR truth/training/known-sites bundle is shared across samples.
    GATK_RESOURCE_BUNDLE()
    resource_bundle = GATK_RESOURCE_BUNDLE.out.resources.first()
    resource_vcfs = resource_bundle.map { omni, omni_tbi, phase, phase_tbi, dbsnp, dbsnp_tbi, hapmap, hapmap_tbi ->
        [hapmap, omni, phase, dbsnp]
    }
    resource_tbis = resource_bundle.map { omni, omni_tbi, phase, phase_tbi, dbsnp, dbsnp_tbi, hapmap, hapmap_tbi ->
        [hapmap_tbi, omni_tbi, phase_tbi, dbsnp_tbi]
    }
    resource_labels = Channel.value([
            '--resource:hapmap,known=false,training=true,truth=true,prior=15.0 hapmap_3.3.hg38.vcf.gz',
            '--resource:omni,known=false,training=true,truth=false,prior=12.0 1000G_omni2.5.hg38.vcf.gz',
            '--resource:1000G,known=false,training=true,truth=false,prior=10.0 1000G_phase1.snps.high_confidence.hg38.vcf.gz',
            '--resource:dbsnp,known=true,training=false,truth=false,prior=2.0 Homo_sapiens_assembly38.dbsnp138.vcf.gz'
    ])

    haplotypecaller_vcf = GATK4_HAPLOTYPECALLER.out.vcf.map { meta, vcf -> tuple(meta.id, meta, vcf) }
    haplotypecaller_tbi = GATK4_HAPLOTYPECALLER.out.tbi.map { meta, tbi -> tuple(meta.id, tbi) }
    // GATK emits the VCF and index on separate channels; match them strictly by
    // sample here because VariantRecalibrator requires both files together.
    recalibrator_input = haplotypecaller_vcf
        .join(haplotypecaller_tbi, failOnMismatch: true, failOnDuplicate: true)
        .map { sample_id, meta, vcf, tbi -> tuple(meta, vcf, tbi) }

    GATK4_VARIANTRECALIBRATOR(
        recalibrator_input,
        resource_vcfs,
        resource_tbis,
        resource_labels,
        reference.map { fasta, fai, dict -> fasta },
        reference.map { fasta, fai, dict -> fai },
        reference.map { fasta, fai, dict -> dict }
    )

    vcf_for_apply = GATK4_HAPLOTYPECALLER.out.vcf.map { meta, vcf -> tuple(meta.id, meta, vcf) }
    tbi_for_apply = GATK4_HAPLOTYPECALLER.out.tbi.map { meta, tbi -> tuple(meta.id, tbi) }
    recal_for_apply = GATK4_VARIANTRECALIBRATOR.out.recal.map { meta, recal -> tuple(meta.id, recal) }
    recal_index_for_apply = GATK4_VARIANTRECALIBRATOR.out.idx.map { meta, idx -> tuple(meta.id, idx) }
    tranches_for_apply = GATK4_VARIANTRECALIBRATOR.out.tranches.map { meta, tranches -> tuple(meta.id, tranches) }

    // ApplyVQSR needs five independently emitted artifacts from the caller and
    // recalibrator. Add each one by sample so an absent or duplicate file fails
    // at the exact channel boundary that supplied it.
    apply_input = vcf_for_apply
        .join(tbi_for_apply, failOnMismatch: true, failOnDuplicate: true)
        .join(recal_for_apply, failOnMismatch: true, failOnDuplicate: true)
        .join(recal_index_for_apply, failOnMismatch: true, failOnDuplicate: true)
        .join(tranches_for_apply, failOnMismatch: true, failOnDuplicate: true)
        .map { sample_id, meta, vcf, tbi, recal, recal_index, tranches ->
            tuple(meta, vcf, tbi, recal, recal_index, tranches)
        }

    GATK4_APPLYVQSR(
        apply_input,
        reference.map { fasta, fai, dict -> fasta },
        reference.map { fasta, fai, dict -> fai },
        reference.map { fasta, fai, dict -> dict }
    )

    results = GATK4_APPLYVQSR.out.vcf.map { meta, vcf -> tuple(meta.id, vcf) }

    emit:
    vcf = results
}
