// Merge SV evidence sources into a final prioritized SV VCF.
process sv_prio {
    container 'beoungl/docker_test:longphase_0.4.0'

    input:
    tuple val(out_prefix), path(sv_phenosv_evidence), path(annovar_sv_vcf), path(hpo_path), val(age_of_onset)

    val(inheritance_mode)
    val(include_clinvar_report)
    val(allow_unphased_comphet)

    output:
    tuple val(out_prefix), path("${out_prefix}.prio.vcf"), emit: prio_vcf
    tuple val(out_prefix), path("${out_prefix}.prio_gene.vcf"), emit: prio_gene_vcf
    tuple val(out_prefix), path("${out_prefix}.frequency_audit.tsv"), emit: frequency_audit

    script:
    def min_score = params.phenosv_score ? params.phenosv_score.toString().trim() : '0.50'
    def gene_filter = params.gene ? params.gene.toString().trim() : ''
    def gene_filter_file = gene_filter ? new File(gene_filter) : null
    if (gene_filter_file?.isFile()) {
        gene_filter = gene_filter_file.readLines().collect { gene -> gene.trim() }.findAll { gene -> gene && !gene.startsWith('#') }.join(',')
    }
    def gnomad_af_ad = params.gnomad_af_ad.toString().trim()
    def gnomad_af_ar = params.gnomad_af_ar.toString().trim()
    def common_sv_af_ad = params.common_sv_af_ad.toString().trim()
    def common_sv_af_ar = params.common_sv_af_ar.toString().trim()
    def common_sv_filter_mode = params.common_sv_filter.toString().trim().toLowerCase()

	"""
	# Convert PhenoSV evidence into inheritance-aware variant- and gene-level reports.

		bash /phenosv_vcf_and_tsv.sh $sv_phenosv_evidence $annovar_sv_vcf $annovar_sv_vcf $out_prefix
	mv ${out_prefix}.phenosv.vcf ${out_prefix}.phenosv.unfiltered.vcf
	python3 /filter_phenosv_vcf.py ${out_prefix}.phenosv.unfiltered.vcf ${out_prefix}.phenosv.vcf --min-score $min_score --genes "$gene_filter"

	python3 /assign_dom_or_rec_sv_only.py ${out_prefix}.phenosv.vcf "$hpo_path" "$age_of_onset" "$inheritance_mode" ${out_prefix}.assigned.vcf --phenosv-score $min_score --genes "$gene_filter"

	python3 /prio_gene_only.py ${out_prefix}.assigned.vcf ${out_prefix}.prio_gene.vcf gene --include-clinvar $include_clinvar_report --allow-unphased-comphet $allow_unphased_comphet --genes "$gene_filter" --gnomad-af-ad ${gnomad_af_ad} --gnomad-af-ar ${gnomad_af_ar} --common-sv-af-ad ${common_sv_af_ad} --common-sv-af-ar ${common_sv_af_ar} --common-sv-filter ${common_sv_filter_mode} --frequency-audit ${out_prefix}.frequency_audit.tsv
	python3 /prio_gene_only.py ${out_prefix}.assigned.vcf ${out_prefix}.prio.vcf variant --include-clinvar $include_clinvar_report --allow-unphased-comphet $allow_unphased_comphet --genes "$gene_filter" --gnomad-af-ad ${gnomad_af_ad} --gnomad-af-ar ${gnomad_af_ar} --common-sv-af-ad ${common_sv_af_ad} --common-sv-af-ar ${common_sv_af_ar} --common-sv-filter ${common_sv_filter_mode}
	"""

}
