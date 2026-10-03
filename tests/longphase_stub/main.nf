nextflow.enable.dsl = 2

include { LONGPHASE_PROCESSING } from '../../subworkflows/longphase_processing'

workflow {
    fixture_root = file("${projectDir}/fixtures")

    records = Channel.of(tuple(
        'sample1',
        file("${fixture_root}/sample.rankvar.tsv"),
        file("${fixture_root}/sample.rankscore.tsv"),
        file("${fixture_root}/sample.clinvar.tsv"),
        file("${fixture_root}/sample.phenosv.tsv"),
        file("${fixture_root}/sample.sv.vcf"),
        file("${fixture_root}/sample.sv_annotation.vcf"),
        file("${fixture_root}/sample.snv.vcf"),
        file("${fixture_root}/sample.bam"),
        file("${fixture_root}/sample.bam.bai"),
        file("${fixture_root}/sample.hpo.txt"),
        'adult'
    ))

    reference_bundle = Channel.value(tuple(
        file("${fixture_root}/reference.fa"),
        file("${fixture_root}/reference.fa.fai")
    ))

    result = LONGPHASE_PROCESSING(
        records,
        reference_bundle,
        Channel.value(false),
        Channel.value('unknown'),
        Channel.value('yes'),
        Channel.value(false)
    )

    checked_prio = result.prio_vcf.map { sample_id, output ->
        assert sample_id == 'sample1'
        assert output.name == 'sample1.prio.vcf'
        "prio:${sample_id}"
    }
    checked_gene = result.prio_gene.map { sample_id, output ->
        assert sample_id == 'sample1'
        assert output.name == 'sample1.prio_gene.vcf'
        "gene:${sample_id}"
    }
    checked_audit = result.frequency_audit.map { sample_id, output ->
        assert sample_id == 'sample1'
        assert output.name == 'sample1.frequency_audit.tsv'
        "audit:${sample_id}"
    }
    checked_haplotag = result.haplotag.map { sample_id, output ->
        assert sample_id == 'sample1'
        assert output.name == 'sample1_haplotag.bam'
        "haplotag:${sample_id}"
    }

    checked_prio
        .mix(checked_gene, checked_audit, checked_haplotag)
        .collect()
        .map { observed ->
            assert observed.sort() == [
                'audit:sample1',
                'gene:sample1',
                'haplotag:sample1',
                'prio:sample1'
            ]
            'LongPhase stub contract OK'
        }
        .view()
}
