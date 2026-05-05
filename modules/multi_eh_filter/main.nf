
// Batch filter of ExpansionHunter outputs for pathogenic repeat expansions.
process multi_eh_filter {
	container ='beoungl/docker_test:eh_filter'


        input:
	tuple val(out_prefix), path(eh_output)

	output:
	path "${out_prefix}.eh.tsv"

	script:

	"""
	python3 /filter_eh/filter_eh.py "$eh_output" > "${out_prefix}.eh.tsv"
	"""


}


