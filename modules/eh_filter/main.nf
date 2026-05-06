
// Filter ExpansionHunter JSON output to loci above pathogenic repeat thresholds.
process eh_filter {
	container ='beoungl/docker_test:eh_filter_0.1'


        input:
	val out_prefix
	path eh_output

	output:
	path "${out_prefix}.eh.tsv"

	script:

	"""
	python3 /filter_eh/filter_eh.py "$eh_output" > "${out_prefix}.eh.tsv"
	"""


}

