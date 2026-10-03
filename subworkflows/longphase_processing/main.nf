include { LONGPHASE_CALL; LONGPHASE_PRIORITIZE } from '../../modules/local/longphase_support'

// Shared keyed LongPhase route. The caller uses the verified historical 1.7.3
// command contract; PipeVar-specific evidence aggregation remains separate.
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
    reference_for_longphase = reference.map { fasta, fai -> tuple([id: 'reference'], fasta, fai) }

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

    call_input = context_by_sample.map { sample_id, context ->
        tuple(context.meta, context.bam, context.bai, context.snv_vcf, context.sv_phase_vcf)
    }
    LONGPHASE_CALL(call_input, reference_for_longphase)

    phased = LONGPHASE_CALL.out.calls.map { meta, phased_snv, phased_sv, haplotag_bam ->
        tuple(meta.id, meta, phased_snv, phased_sv, haplotag_bam)
    }

    priority_input = context_by_sample
        .join(phased, failOnMismatch: true, failOnDuplicate: true)

    unified_priority_input = priority_input.map {
        _sample_id, context, _call_meta, phased_snv, phased_sv, _haplotag_bam ->
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

    haplotag = LONGPHASE_CALL.out.calls.map { meta, _phased_snv, _phased_sv, bam -> tuple(meta.id, bam) }

    emit:
    prio_vcf
    prio_gene
    frequency_audit
    haplotag
}
