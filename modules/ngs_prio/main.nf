
// Integrate NGS SNV/SV evidence and produce a prioritized VCF.
process ngs_prio {
	container ='beoungl/docker_test:longphase_0.2.31'

	input:
	val(out_prefix)
	path(snv_rankvar)
	path(snv_rankscore)
	path(snv_pathogenic)
	path(sv_pathogenic)
	path(sv_vcf_path)
	path(snv_vcf_path)
	path(hpo_path)
	val(inheritance_mode)
	val(include_clinvar_report)
	val(allow_unphased_comphet)

	output:
	path "${out_prefix}.prio.vcf"
	path "${out_prefix}.prio_gene.vcf"

	
	script:
	def min_score = params.phenosv_score ?: '0.50'
	def gene_filter = params.gene ?: ''
	def sv_only = params.prioritize_sv_only ?: 'no'


	"""

	#No bam file processing, so only sv + snp processing.

	        bash /phenosv_vcf_and_tsv.sh $sv_pathogenic $sv_vcf_path $sv_vcf_path $out_prefix
	mv ${out_prefix}.phenosv.vcf ${out_prefix}.phenosv.unfiltered.vcf
	python3 /filter_phenosv_vcf.py ${out_prefix}.phenosv.unfiltered.vcf ${out_prefix}.phenosv.vcf --min-score $min_score --genes "$gene_filter"

	if [[ "$sv_only" == "yes" ]]; then
	    python3 /assign_dom_or_rec_sv_only.py ${out_prefix}.phenosv.vcf $hpo_path $inheritance_mode ${out_prefix}.assigned.vcf --phenosv-score $min_score --genes "$gene_filter"
	else
	    bash /clinvar_vcf_and_txt.sh $snv_pathogenic $snv_vcf_path $out_prefix
	    bash /rankscore_vcf_and_txt.sh $snv_rankscore $snv_vcf_path $out_prefix
	    bash /rankvar_vcf_and_tsv.sh $snv_rankvar $snv_vcf_path $out_prefix
	    python3 /assign_dom_or_rec.py ${out_prefix}.clinvar.vcf ${out_prefix}.phenosv.vcf ${out_prefix}.rankscore.vcf ${out_prefix}.rankvar.vcf $hpo_path $inheritance_mode ${out_prefix}.assigned.vcf --phenosv-score $min_score --genes "$gene_filter"
	fi

	# Keep the ranked gene-level report and also emit a ranked variant-level VCF.
	python3 /prio_gene_only.py ${out_prefix}.assigned.vcf ${out_prefix}.prio_gene.vcf gene --include-clinvar $include_clinvar_report --allow-unphased-comphet $allow_unphased_comphet --genes "$gene_filter"
	python3 /prio_gene_only.py ${out_prefix}.assigned.vcf ${out_prefix}.prio.vcf variant --include-clinvar $include_clinvar_report --allow-unphased-comphet $allow_unphased_comphet --genes "$gene_filter"

	"""

}
