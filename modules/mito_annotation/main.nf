process mito_annotation {
	container = 'beoungl/docker_test:mito_annotation_0.1'

	input:
	path mito_vcf
	val out_prefix

	output:
	path "${out_prefix}.mito.annotated.tsv"
	path "${out_prefix}.mito.annotated.vcf.gz"
	path "${out_prefix}.mito.annotated.vcf.gz.tbi"

	script:
	"""
	set -euo pipefail

	python3 /opt/mito/bin/annotate_mito_variants.py \
	    --vcf "$mito_vcf" \
	    --out-prefix "$out_prefix" \
	    --hmtvar-status "${params.hmtvar_data ? 'bundled' : 'asset_missing'}"
	"""
}
