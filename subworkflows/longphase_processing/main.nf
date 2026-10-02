include { LONGPHASE_PREPARE_INPUT; LONGPHASE_PRIORITIZE } from '../../modules/local/longphase_support'
include { LONGPHASE_PHASE } from '../../modules/nf-core/longphase/phase'
include { LONGPHASE_HAPLOTAG } from '../../modules/nf-core/longphase/haplotag'

// Shared keyed LongPhase route. Generic phase/haplotag operations use pinned
// nf-core modules; PipeVar-specific evidence aggregation remains local.
workflow LONGPHASE_PROCESSING {
    take:
    records
    reference_bundle
    batch_input
    inheritance_mode
    include_clinvar_report
    allow_unphased_comphet

    main:
    reference = reference_bundle
    reference_fasta = reference.map { fasta, fai -> tuple([id: 'reference'], fasta) }
    reference_fai = reference.map { fasta, fai -> tuple([id: 'reference'], fai) }

    // The incoming evidence record is wide to preserve the established process tuple contract.
    // Convert it once to a named context so later joins remain reviewable.
    context_by_sample = records.map {
        sample_id,
        snv_rankvar,
        snv_rankscore,
        snv_clinvar_evidence,
        sv_phenosv_evidence,
        sv_phase_vcf,
        sv_annotation_vcf,
        snv_vcf,
        bam,
        bai,
        hpo,
        age ->
        tuple(
            sample_id,
            [
                meta: [id: sample_id, single_end: false, batch_input: batch_input],
                snv_rankvar: snv_rankvar,
                snv_rankscore: snv_rankscore,
                snv_clinvar_evidence: snv_clinvar_evidence,
                sv_phenosv_evidence: sv_phenosv_evidence,
                sv_phase_vcf: sv_phase_vcf,
                sv_annotation_vcf: sv_annotation_vcf,
                snv_vcf: snv_vcf,
                bam: bam,
                bai: bai,
                hpo: hpo,
                age: age
            ]
        )
    }

    prepare_input = context_by_sample.map { sample_id, context ->
        tuple(context.meta, context.snv_vcf)
    }
    LONGPHASE_PREPARE_INPUT(
        prepare_input,
        reference.map { fasta, fai -> fai }
    )

    prepared_by_id = LONGPHASE_PREPARE_INPUT.out.vcf.map { meta, prepared_vcf ->
        tuple(meta.id, prepared_vcf)
    }
    // (meta, BAM, BAI, small-variant VCF, SV VCF, regions); [] phases the whole reference.
    phase_input = context_by_sample
        .join(prepared_by_id, failOnMismatch: true, failOnDuplicate: true)
        .map { sample_id, context, prepared_vcf ->
            tuple(context.meta, context.bam, context.bai, prepared_vcf, context.sv_phase_vcf, [])
        }

    LONGPHASE_PHASE(phase_input, reference_fasta, reference_fai)

    // The nf-core phase process emits phased SNV and SV files separately. Match
    // them by sample once, then reuse the paired result for haplotagging and ranking.
    phased = LONGPHASE_PHASE.out.snv_vcf
        .map { meta, phased_snv -> tuple(meta.id, meta, phased_snv) }
        .join(
            LONGPHASE_PHASE.out.sv_vcf.map { meta, phased_sv -> tuple(meta.id, phased_sv) },
            failOnMismatch: true,
            failOnDuplicate: true
        )

    // Haplotagging and prioritization require different projections of the same
    // independently phased result, so both join back to the retained sample context.
    haplotag_input = context_by_sample
        .map { sample_id, context ->
            tuple(sample_id, context.meta, context.bam, context.bai)
        }
        .join(phased, failOnMismatch: true, failOnDuplicate: true)
        .map { sample_id, meta, bam, bai, _phase_meta, phased_snv, phased_sv ->
            tuple(meta, bam, bai, phased_snv, phased_sv, [])
        }
    LONGPHASE_HAPLOTAG(haplotag_input, reference_fasta, reference_fai)

    priority_input = context_by_sample
        .join(phased, failOnMismatch: true, failOnDuplicate: true)

    unified_priority_input = priority_input.map {
        _sample_id, context, _phase_meta, phased_snv, phased_sv ->
        tuple(
            context.meta,
            phased_snv,
            phased_sv,
            context.sv_annotation_vcf,
            context.sv_phenosv_evidence,
            context.snv_rankscore,
            context.snv_clinvar_evidence,
            context.snv_rankvar,
            context.hpo,
            context.age
        )
    }
    LONGPHASE_PRIORITIZE(
        unified_priority_input,
        inheritance_mode,
        include_clinvar_report,
        allow_unphased_comphet
    )
    prio_vcf = LONGPHASE_PRIORITIZE.out.prio_vcf.map { meta, path -> tuple(meta.id, path) }
    prio_gene = LONGPHASE_PRIORITIZE.out.prio_gene_vcf.map { meta, path -> tuple(meta.id, path) }
    frequency_audit = LONGPHASE_PRIORITIZE.out.frequency_audit.map { meta, path -> tuple(meta.id, path) }

    haplotag = LONGPHASE_HAPLOTAG.out.bam.map { meta, bam -> tuple(meta.id, bam) }

    emit:
    prio_vcf
    prio_gene
    frequency_audit
    haplotag
}
