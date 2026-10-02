// Reheader and validate the SNV VCF before handing it to the vendored LongPhase module.
process LONGPHASE_PREPARE_INPUT {
    tag "${meta.id}"
    label 'process_low'

    container 'community.wave.seqera.io/library/htslib_samtools:1.23.1--5b6bb4ede7e612e5'

    input:
    tuple val(meta), path(snv_vcf)
    path fai

    output:
    tuple val(meta), path("${meta.id}.longphase_input.vcf"), emit: vcf

    script:
    """
    bcftools reheader \\
        --fai ${fai} \\
        --output ${meta.id}.longphase_input.vcf \\
        ${snv_vcf}

    awk '
    BEGIN {
        while ((getline line < ARGV[1]) > 0) {
            split(line, columns, "\t")
            reference_length[columns[1]] = columns[2]
        }
        close(ARGV[1])
        ARGV[1] = ""
    }
    /^##contig=</ {
        contig_headers++
        if (\$0 !~ /,length=[1-9][0-9]*([,>])/) {
            print "ERROR: LongPhase SNV VCF contig header lacks a positive numeric length: " \$0 > "/dev/stderr"
            errors = 1
        }
        next
    }
    /^#/ { next }
    {
        if (!(\$1 in reference_length)) {
            print "ERROR: LongPhase SNV VCF contains contig absent from reference FAI at " \$1 ":" \$2 > "/dev/stderr"
            errors = 1
        }
        if (\$2 !~ /^[1-9][0-9]*\$/) {
            print "ERROR: LongPhase SNV VCF contains an invalid POS at " \$1 ":" \$2 > "/dev/stderr"
            errors = 1
        }
        else if ((\$1 in reference_length) && \$2 > reference_length[\$1]) {
            print "ERROR: LongPhase SNV VCF POS exceeds reference contig length at " \$1 ":" \$2 > "/dev/stderr"
            errors = 1
        }
    }
    END {
        if (contig_headers == 0) {
            print "ERROR: LongPhase SNV VCF has no contig headers after FAI reheadering" > "/dev/stderr"
            errors = 1
        }
        exit errors
    }' ${fai} ${meta.id}.longphase_input.vcf
    """

    stub:
    """
    touch ${meta.id}.longphase_input.vcf
    """
}

// Prioritize LongPhase evidence with optional age-aware inheritance scoring.
process LONGPHASE_PRIORITIZE {
    tag "${meta.id}"
    label 'process_low'

    container 'beoungl/docker_test:longphase_0.4.0'

    input:
    tuple val(meta), path(phased_snv), path(phased_sv), path(sv_annotation_vcf), path(sv_phenosv_evidence), path(snv_rankscore), path(snv_clinvar_evidence), path(snv_rankvar), path(hpo_path), val(age_of_onset)
    val inheritance_mode
    val include_clinvar_report
    val allow_unphased_comphet

    output:
    tuple val(meta), path("${meta.id}.prio.vcf"), emit: prio_vcf
    tuple val(meta), path("${meta.id}.prio_gene.vcf"), emit: prio_gene_vcf
    tuple val(meta), path("${meta.id}.frequency_audit.tsv"), emit: frequency_audit

    script:
    def min_score = params.phenosv_score ?: '0.50'
    def gene_filter = params.gene ?: ''
    def sv_only = params.prioritize_sv_only ?: 'no'
    """
    bcftools view -Ov -o ${meta.id}_phased.vcf ${phased_snv}
    bcftools view -Ov -o ${meta.id}_phased_SV.vcf ${phased_sv}

    bash /phenosv_vcf_and_tsv.sh ${sv_phenosv_evidence} ${meta.id}_phased_SV.vcf ${sv_annotation_vcf} ${meta.id}
    mv ${meta.id}.phenosv.vcf ${meta.id}.phenosv.unfiltered.vcf
    python3 /filter_phenosv_vcf.py ${meta.id}.phenosv.unfiltered.vcf ${meta.id}.phenosv.vcf --min-score ${min_score} --genes "${gene_filter}"

    if [[ "${sv_only}" == "yes" ]]; then
        python3 /assign_dom_or_rec_sv_only.py ${meta.id}.phenosv.vcf "${hpo_path}" "${age_of_onset}" "${inheritance_mode}" ${meta.id}.assigned.vcf --phenosv-score ${min_score} --genes "${gene_filter}"
    else
        bash /clinvar_vcf_and_txt.sh ${snv_clinvar_evidence} ${meta.id}_phased.vcf ${meta.id}
        bash /rankscore_vcf_and_txt.sh ${snv_rankscore} ${meta.id}_phased.vcf ${meta.id}
        bash /rankvar_vcf_and_tsv.sh ${snv_rankvar} ${meta.id}_phased.vcf ${meta.id}
        python3 /assign_dom_or_rec.py ${meta.id}.clinvar.vcf ${meta.id}.phenosv.vcf ${meta.id}.rankscore.vcf ${meta.id}.rankvar.vcf "${hpo_path}" "${age_of_onset}" "${inheritance_mode}" ${meta.id}.assigned.vcf --phenosv-score ${min_score} --genes "${gene_filter}"
    fi

    python3 /prio_gene_only.py ${meta.id}.assigned.vcf ${meta.id}.prio_gene.vcf gene --include-clinvar ${include_clinvar_report} --allow-unphased-comphet ${allow_unphased_comphet} --genes "${gene_filter}" --gnomad-af-ad ${params.gnomad_af_ad} --gnomad-af-ar ${params.gnomad_af_ar} --common-sv-af-ad ${params.common_sv_af_ad} --common-sv-af-ar ${params.common_sv_af_ar} --common-sv-filter ${params.common_sv_filter} --frequency-audit ${meta.id}.frequency_audit.tsv
    python3 /prio_gene_only.py ${meta.id}.assigned.vcf ${meta.id}.prio.vcf variant --include-clinvar ${include_clinvar_report} --allow-unphased-comphet ${allow_unphased_comphet} --genes "${gene_filter}" --gnomad-af-ad ${params.gnomad_af_ad} --gnomad-af-ar ${params.gnomad_af_ar} --common-sv-af-ad ${params.common_sv_af_ad} --common-sv-af-ar ${params.common_sv_af_ar} --common-sv-filter ${params.common_sv_filter}
    """

    stub:
    """
    touch ${meta.id}.prio.vcf
    touch ${meta.id}.prio_gene.vcf
    touch ${meta.id}.frequency_audit.tsv
    """
}
