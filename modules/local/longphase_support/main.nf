// Run the verified historical LongPhase 1.7.3 command contract behind the
// simplified metadata-keyed interface. Phasing and haplotagging intentionally
// stay in one image so their binary and utility requirements cannot drift.
process LONGPHASE_CALL {
    tag "${meta.id}"
    label 'process_high'

    container 'beoungl/docker_test:longphase_0.4.0'

    input:
    tuple val(meta), path(bam), path(bai), path(snv_vcf), path(sv_vcf)
    tuple val(reference_meta), path(fasta), path(fai)

    output:
    tuple val(meta), path("${meta.id}_phased.vcf"), path("${meta.id}_phased_SV.vcf"), path("${meta.id}_haplotag.bam"), emit: calls

    script:
    def platform_args = task.ext.args != null ? task.ext.args : '--ont'
    def phase_threads = task.ext.phase_threads != null ? task.ext.phase_threads : task.cpus
    def haplotag_threads = task.ext.haplotag_threads != null ? task.ext.haplotag_threads : task.cpus
    """
    echo "PIPEVAR_LONGPHASE_RESOURCES process=${task.process} attempt=${task.attempt} cpus=${task.cpus} memory='${task.memory}' phase_threads=${phase_threads} haplotag_threads=${haplotag_threads}" >&2
    longphase_version=\$(/longphase_linux-x64 --version 2>&1 | sed -n '1p')
    echo "PIPEVAR_LONGPHASE_VERSION \${longphase_version}" >&2
    if [[ -r /sys/fs/cgroup/memory.max ]]; then
        echo "PIPEVAR_LONGPHASE_CGROUP memory.max=\$(cat /sys/fs/cgroup/memory.max)" >&2
    elif [[ -r /sys/fs/cgroup/memory/memory.limit_in_bytes ]]; then
        echo "PIPEVAR_LONGPHASE_CGROUP memory.limit_in_bytes=\$(cat /sys/fs/cgroup/memory/memory.limit_in_bytes)" >&2
    else
        echo "PIPEVAR_LONGPHASE_CGROUP unavailable" >&2
    fi

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

    /longphase_linux-x64 phase \\
        -s ${meta.id}.longphase_input.vcf \\
        --sv-file=${sv_vcf} \\
        -t ${phase_threads} \\
        -o ${meta.id}_phased \\
        ${platform_args} \\
        -b ${bam} \\
        -r ${fasta}

    /longphase_linux-x64 haplotag \\
        -r ${fasta} \\
        -s ${meta.id}_phased.vcf \\
        --sv-file ${meta.id}_phased_SV.vcf \\
        -b ${bam} \\
        -t ${haplotag_threads} \\
        -o ${meta.id}_haplotag
    """

    stub:
    """
    touch ${meta.id}_phased.vcf
    touch ${meta.id}_phased_SV.vcf
    touch ${meta.id}_haplotag.bam
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
    def min_score = params.phenosv_score ? params.phenosv_score.toString().trim() : '0.50'
    def gene_filter = params.gene ? params.gene.toString().trim() : ''
    def gene_filter_file = gene_filter ? new File(gene_filter) : null
    if (gene_filter_file?.isFile()) {
        gene_filter = gene_filter_file.readLines().collect { gene -> gene.trim() }.findAll { gene -> gene && !gene.startsWith('#') }.join(',')
    }
    def sv_only = params.prioritize_sv_only ? params.prioritize_sv_only.toString().trim().toLowerCase() : 'no'
    def gnomad_af_ad = params.gnomad_af_ad.toString().trim()
    def gnomad_af_ar = params.gnomad_af_ar.toString().trim()
    def common_sv_af_ad = params.common_sv_af_ad.toString().trim()
    def common_sv_af_ar = params.common_sv_af_ar.toString().trim()
    def common_sv_filter_mode = params.common_sv_filter.toString().trim().toLowerCase()
    """
    bcftools view -Ov -o ${meta.id}.priority_input.snv.vcf ${phased_snv}
    bcftools view -Ov -o ${meta.id}.priority_input.sv.vcf ${phased_sv}

    bash /phenosv_vcf_and_tsv.sh ${sv_phenosv_evidence} ${meta.id}.priority_input.sv.vcf ${sv_annotation_vcf} ${meta.id}
    mv ${meta.id}.phenosv.vcf ${meta.id}.phenosv.unfiltered.vcf
    python3 /filter_phenosv_vcf.py ${meta.id}.phenosv.unfiltered.vcf ${meta.id}.phenosv.vcf --min-score ${min_score} --genes "${gene_filter}"

    if [[ "${sv_only}" == "yes" ]]; then
        python3 /assign_dom_or_rec_sv_only.py ${meta.id}.phenosv.vcf "${hpo_path}" "${age_of_onset}" "${inheritance_mode}" ${meta.id}.assigned.vcf --phenosv-score ${min_score} --genes "${gene_filter}"
    else
        bash /clinvar_vcf_and_txt.sh ${snv_clinvar_evidence} ${meta.id}.priority_input.snv.vcf ${meta.id}
        bash /rankscore_vcf_and_txt.sh ${snv_rankscore} ${meta.id}.priority_input.snv.vcf ${meta.id}
        bash /rankvar_vcf_and_tsv.sh ${snv_rankvar} ${meta.id}.priority_input.snv.vcf ${meta.id}
        python3 /assign_dom_or_rec.py ${meta.id}.clinvar.vcf ${meta.id}.phenosv.vcf ${meta.id}.rankscore.vcf ${meta.id}.rankvar.vcf "${hpo_path}" "${age_of_onset}" "${inheritance_mode}" ${meta.id}.assigned.vcf --phenosv-score ${min_score} --genes "${gene_filter}"
    fi

    python3 /prio_gene_only.py ${meta.id}.assigned.vcf ${meta.id}.prio_gene.vcf gene --include-clinvar ${include_clinvar_report} --allow-unphased-comphet ${allow_unphased_comphet} --genes "${gene_filter}" --gnomad-af-ad ${gnomad_af_ad} --gnomad-af-ar ${gnomad_af_ar} --common-sv-af-ad ${common_sv_af_ad} --common-sv-af-ar ${common_sv_af_ar} --common-sv-filter ${common_sv_filter_mode} --frequency-audit ${meta.id}.frequency_audit.tsv
    python3 /prio_gene_only.py ${meta.id}.assigned.vcf ${meta.id}.prio.vcf variant --include-clinvar ${include_clinvar_report} --allow-unphased-comphet ${allow_unphased_comphet} --genes "${gene_filter}" --gnomad-af-ad ${gnomad_af_ad} --gnomad-af-ar ${gnomad_af_ar} --common-sv-af-ad ${common_sv_af_ad} --common-sv-af-ar ${common_sv_af_ar} --common-sv-filter ${common_sv_filter_mode}
    """

    stub:
    """
    touch ${meta.id}.prio.vcf
    touch ${meta.id}.prio_gene.vcf
    touch ${meta.id}.frequency_audit.tsv
    """
}
