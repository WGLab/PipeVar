process multi_mito_clair3_postprocess {
	container = 'beoungl/docker_test:mito_clair3_postprocess_0.1'

	input:
	tuple val(out_prefix), path(raw_vcf)

	output:
	tuple val(out_prefix), path("${out_prefix}.mito.vcf.gz"), path("${out_prefix}.mito.vcf.gz.tbi")

	script:
	"""
	bcftools norm -m -any "$raw_vcf" -Ov -o "${out_prefix}.mito.split.vcf"

	python3 /usr/local/bin/clair3_mito_adapt.py \\
	    --input "${out_prefix}.mito.split.vcf" \\
	    --output "${out_prefix}.mito.raw.vcf"

	bgzip -f -c "${out_prefix}.mito.raw.vcf" > "${out_prefix}.mito.vcf.gz"
	tabix -f -p vcf "${out_prefix}.mito.vcf.gz"
	rm -f "${out_prefix}.mito.split.vcf" "${out_prefix}.mito.raw.vcf"
	"""
}
