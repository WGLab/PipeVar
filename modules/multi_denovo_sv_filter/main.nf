// Batch de novo filtering for keyed ANNOVAR structural-variant VCF outputs.
process multi_denovo_sv_filter {
	container = 'beoungl/docker_test:denovo_snv_sv_filter_0.2'
	tag 'sv'

	input:
	tuple val(out_prefixes), path(annovar_sv_vcfs)
	path pedigree_csv
	val role_column
	val family_column
	val vcf_sample_column
	val exclude_contigs
	val sv_min_reciprocal_overlap

	output:
	path "*.denovo.sv.exonic.vcf", emit: vcf
	path "denovo.sv.summary.tsv", emit: summary

	script:
	def prefixes = out_prefixes instanceof List ? out_prefixes : [out_prefixes]
	def vcf_list = annovar_sv_vcfs instanceof List ? annovar_sv_vcfs : [annovar_sv_vcfs]
	if (prefixes.size() != vcf_list.size()) {
		throw new IllegalArgumentException("multi_denovo_sv_filter requires matching out_prefix and VCF counts")
	}
	def input_args = (0..<prefixes.size()).collect { idx ->
		"--sv-input '${prefixes[idx]}' '${vcf_list[idx]}'"
	}.join(' ')
	"""
	/usr/local/bin/filter_denovo_snv_sv.py \\
	    --mode sv \\
	    --pedigree-csv $pedigree_csv \\
	    --role-column '$role_column' \\
	    --family-column '$family_column' \\
	    --vcf-sample-column '$vcf_sample_column' \\
	    --exclude-contigs '$exclude_contigs' \\
	    --sv-min-reciprocal-overlap '$sv_min_reciprocal_overlap' \\
	    $input_args
	"""
}
