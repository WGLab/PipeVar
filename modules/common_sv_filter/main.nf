// Filter common structural variants from a single ANNOVAR SV VCF.
process common_sv_filter {
	container = 'beoungl/docker_test:common_sv_filter_0.4'

	input:
	path vcf
	val out_prefix

	output:
	path "${out_prefix}.common_sv_filtered.vcf", emit: filtered_vcf
	path "${out_prefix}.common_sv_removed.vcf", emit: removed_vcf
	path "${out_prefix}.common_sv_filter.summary.tsv", emit: summary

	script:
	def common_af = params.common_sv_af ?: '0.01'
	def reciprocal_overlap = params.common_sv_reciprocal_overlap ?: '0.5'
	def distance = params.common_sv_distance ?: '1000'
	def ins_distance = params.common_sv_ins_distance ?: '500'
	def ins_identity = params.common_sv_ins_identity ?: '0.8'
	"""
	python3 /usr/local/bin/filter_common_svs.py \\
	    --input-vcf $vcf \\
	    --output-vcf ${out_prefix}.common_sv_filtered.vcf \\
	    --removed-vcf ${out_prefix}.common_sv_removed.vcf \\
	    --summary-tsv ${out_prefix}.common_sv_filter.summary.tsv \\
	    --common-af "$common_af" \\
	    --reciprocal-overlap "$reciprocal_overlap" \\
	    --distance "$distance" \\
	    --ins-distance "$ins_distance" \\
	    --ins-identity "$ins_identity"
	"""
}
