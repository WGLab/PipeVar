// Batch de novo filtering for keyed raw/called small-variant VCF outputs.
process multi_denovo_snv_vcf_filter {
	container = 'beoungl/docker_test:denovo_snv_sv_filter_0.3'
	tag 'snv-vcf'

	input:
	tuple val(out_prefixes), path(snv_vcfs)
	path pedigree_csv
	val role_column
	val family_column
	val vcf_sample_column
	val exclude_contigs

	output:
	path "*.denovo.snv.vcf", emit: vcf
	path "denovo.snv.summary.tsv", emit: summary
	path "denovo.snv.bindings.tsv", emit: bindings

	script:
	def prefixes = out_prefixes instanceof List ? out_prefixes : [out_prefixes]
	def vcf_list = snv_vcfs instanceof List ? snv_vcfs : [snv_vcfs]
	if (prefixes.size() != vcf_list.size()) {
		throw new IllegalArgumentException("multi_denovo_snv_vcf_filter requires matching out_prefix and VCF counts")
	}
	def input_args = (0..<prefixes.size()).collect { idx ->
		"--snv-vcf-input '${prefixes[idx]}' '${vcf_list[idx]}'"
	}.join(' ')
	"""
	/usr/local/bin/filter_denovo_snv_sv.py \\
	    --mode snv \\
	    --pedigree-csv $pedigree_csv \\
	    --role-column '$role_column' \\
	    --family-column '$family_column' \\
	    --vcf-sample-column '$vcf_sample_column' \\
	    --exclude-contigs '$exclude_contigs' \\
	    $input_args
	"""
}
