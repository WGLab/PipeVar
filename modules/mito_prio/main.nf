process mito_prio {
	container = 'beoungl/docker_test:mito_annotation_0.4.1'

	input:
	path annotated_tsv
	val out_prefix

	output:
	path "${out_prefix}.mito.prioritized.tsv"

	script:
	"""
	python3 /opt/mito/bin/prioritize_mito_variants.py \
	    --input "$annotated_tsv" \
	    --out "${out_prefix}.mito.prioritized.tsv" \
	    --min-vaf "${params.mito_min_vaf}" \
	    --min-depth "${params.mito_min_depth}" \
	    --min-alt-reads "${params.mito_min_alt_reads}" \
	    --min-priority-vaf "${params.mito_gui_min_af}" \
	    --min-apogee2-score "${params.mito_gui_min_apogee2}" \
	    --min-mitotip-score "${params.mito_gui_min_mitotip}"
	"""
}
