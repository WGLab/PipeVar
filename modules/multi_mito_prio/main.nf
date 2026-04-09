process multi_mito_prio {
	container = 'beoungl/docker_test:mito_annotation_0.1'

	input:
	tuple val(out_prefix), path(annotated_tsv), path(annotated_vcf), path(annotated_tbi)

	output:
	tuple val(out_prefix), path("${out_prefix}.mito.prioritized.tsv")

	script:
	"""
	set -euo pipefail

	python3 /opt/mito/bin/prioritize_mito_variants.py \
	    --input "$annotated_tsv" \
	    --out "${out_prefix}.mito.prioritized.tsv" \
	    --min-vaf "${params.mito_min_vaf}" \
	    --min-depth "${params.mito_min_depth}" \
	    --min-alt-reads "${params.mito_min_alt_reads}"
	"""
}
