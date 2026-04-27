
// Filter ExpansionHunter JSON output to loci above pathogenic repeat thresholds.
process eh_filter {
	container ='beoungl/docker_test:eh_filter'


        input:
	val out_prefix
	path eh_output

	output:
	path "${out_prefix}.eh.tsv"

	script:
	def filtered_input = eh_output.name.endsWith('.gz') ? "${out_prefix}.eh_input.json" : eh_output

	"""
	if [[ "$eh_output" == *.gz ]]; then
		python3 - "$eh_output" "$filtered_input" <<'PY'
import gzip
import sys

with gzip.open(sys.argv[1], "rt") as src, open(sys.argv[2], "w") as dst:
    dst.write(src.read())
PY
	fi
	python3 /filter_eh/filter_eh.py $filtered_input > ${out_prefix}.eh.tsv
	
	
	"""


}


