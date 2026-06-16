// Batch de novo filtering for keyed ANNOVAR small-variant TXT/VCF outputs.
process multi_denovo_snv_filter {
	container = 'beoungl/docker_test:denovo_snv_sv_filter_0.1'
	tag 'snv'

	input:
	tuple val(out_prefixes), path(annovar_txts), path(annovar_vcfs)
	path pedigree_csv
	val role_column
	val family_column
	val vcf_sample_column
	val exclude_contigs

	output:
	path "*.denovo.snv.hg38_multianno.txt", emit: txt
	path "*.denovo.snv.hg38_multianno.vcf", emit: vcf
	path "denovo.snv.summary.tsv", emit: summary

	script:
	def prefixes = out_prefixes instanceof List ? out_prefixes : [out_prefixes]
	def txt_list = annovar_txts instanceof List ? annovar_txts : [annovar_txts]
	def vcf_list = annovar_vcfs instanceof List ? annovar_vcfs : [annovar_vcfs]
	if (prefixes.size() != txt_list.size() || prefixes.size() != vcf_list.size()) {
		throw new IllegalArgumentException("multi_denovo_snv_filter requires matching out_prefix, TXT, and VCF counts")
	}
	def input_args = (0..<prefixes.size()).collect { idx ->
		"--snv-input '${prefixes[idx]}' '${txt_list[idx]}' '${vcf_list[idx]}'"
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
