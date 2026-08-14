
// Merge SV evidence sources into a final prioritized SV VCF.
process sv_prio {
	container ='beoungl/docker_test:longphase_0.4.0'

	input:
	val(out_prefix) 
	path(sv_pathogenic)
	path(annovar_sv_vcf)
	path(hpo_path)

	val(inheritance_mode)
	val(include_clinvar_report)
	val(allow_unphased_comphet)

	output:
	path "${out_prefix}.prio.vcf"
	path "${out_prefix}.prio_gene.vcf"
	path "${out_prefix}.frequency_audit.tsv"
	
	script:
	def min_score = params.phenosv_score ?: '0.50'
	def gene_filter = params.gene ?: ''


	"""
	#Final process to organize the results and show them in ACMG guideline format; might take long to develop.

		bash /phenosv_vcf_and_tsv.sh $sv_pathogenic $annovar_sv_vcf $annovar_sv_vcf $out_prefix
	mv ${out_prefix}.phenosv.vcf ${out_prefix}.phenosv.unfiltered.vcf
	python3 /filter_phenosv_vcf.py ${out_prefix}.phenosv.unfiltered.vcf ${out_prefix}.phenosv.vcf --min-score $min_score --genes "$gene_filter"


	#Get SNP from RankScore + RankVar + ClinVar, and based on score + homozygosity, rank them.

	python3 /assign_dom_or_rec_sv_only.py ${out_prefix}.phenosv.vcf $hpo_path $inheritance_mode ${out_prefix}.assigned.vcf --phenosv-score $min_score --genes "$gene_filter"

	python3 /prio_gene_only.py ${out_prefix}.assigned.vcf ${out_prefix}.prio_gene.vcf gene --include-clinvar $include_clinvar_report --allow-unphased-comphet $allow_unphased_comphet --genes "$gene_filter" --gnomad-af-ad ${params.gnomad_af_ad} --gnomad-af-ar ${params.gnomad_af_ar} --common-sv-af-ad ${params.common_sv_af_ad} --common-sv-af-ar ${params.common_sv_af_ar} --common-sv-filter ${params.common_sv_filter} --de-novo ${params.denovo_filter} --frequency-audit ${out_prefix}.frequency_audit.tsv
	python3 /prio_gene_only.py ${out_prefix}.assigned.vcf ${out_prefix}.prio.vcf variant --include-clinvar $include_clinvar_report --allow-unphased-comphet $allow_unphased_comphet --genes "$gene_filter" --gnomad-af-ad ${params.gnomad_af_ad} --gnomad-af-ar ${params.gnomad_af_ar} --common-sv-af-ad ${params.common_sv_af_ad} --common-sv-af-ar ${params.common_sv_af_ar} --common-sv-filter ${params.common_sv_filter} --de-novo ${params.denovo_filter}
	

	#Potential homozygotes are ranked higher than single heterozygotes, but nothing conclusive.
	#Also gene-based.

	#Use similar strategy as the all long-phase analysis as well here.


	"""

}
