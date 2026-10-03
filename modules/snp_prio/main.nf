// Merge SNV/indel evidence sources into final prioritized small-variant VCFs.
process small_variant_prio {
    container 'beoungl/docker_test:longphase_0.4.0'

    input:
    tuple val(out_prefix), path(snv_rankscore), path(snv_clinvar_evidence), path(snv_rankvar), path(annovar_vcf), path(hpo_path), val(age_of_onset)

    val(inheritance_mode)
    val(include_clinvar_report)
    val(allow_unphased_comphet)

    output:
    tuple val(out_prefix), path("${out_prefix}.prio.vcf"), emit: prio_vcf
    tuple val(out_prefix), path("${out_prefix}.prio_gene.vcf"), emit: prio_gene_vcf
    tuple val(out_prefix), path("${out_prefix}.frequency_audit.tsv"), emit: frequency_audit

    script:

	"""
	# Combine ClinVar, RankScore, and RankVar evidence before inheritance-aware prioritization.

	bash /clinvar_vcf_and_txt.sh $snv_clinvar_evidence $annovar_vcf $out_prefix

        bash /rankscore_vcf_and_txt.sh $snv_rankscore $annovar_vcf $out_prefix

        bash /rankvar_vcf_and_tsv.sh $snv_rankvar $annovar_vcf $out_prefix


	python3 /assign_dom_or_rec_snp_only.py ${out_prefix}.clinvar.vcf ${out_prefix}.rankscore.vcf ${out_prefix}.rankvar.vcf "$hpo_path" "$age_of_onset" "$inheritance_mode" ${out_prefix}.assigned.vcf

	python3 /prio_gene_only.py ${out_prefix}.assigned.vcf ${out_prefix}.prio_gene.vcf gene --include-clinvar $include_clinvar_report --allow-unphased-comphet $allow_unphased_comphet --gnomad-af-ad ${params.gnomad_af_ad} --gnomad-af-ar ${params.gnomad_af_ar} --common-sv-af-ad ${params.common_sv_af_ad} --common-sv-af-ar ${params.common_sv_af_ar} --common-sv-filter ${params.common_sv_filter} --frequency-audit ${out_prefix}.frequency_audit.tsv
	python3 /prio_gene_only.py ${out_prefix}.assigned.vcf ${out_prefix}.prio.vcf variant --include-clinvar $include_clinvar_report --allow-unphased-comphet $allow_unphased_comphet --gnomad-af-ad ${params.gnomad_af_ad} --gnomad-af-ar ${params.gnomad_af_ar} --common-sv-af-ad ${params.common_sv_af_ad} --common-sv-af-ar ${params.common_sv_af_ar} --common-sv-filter ${params.common_sv_filter}
	"""

}
