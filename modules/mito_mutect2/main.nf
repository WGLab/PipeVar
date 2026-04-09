process mito_mutect2 {
	container = 'beoungl/docker_test:mito_mutect2_0.1'

	input:
	tuple path(bam), path(index_file), path(metrics_file)
	val out_prefix
	tuple path(ref_fa), path(fa_index), path(dict_index)
	val mito_contig

	output:
	path "${out_prefix}.mito.vcf.gz"
	path "${out_prefix}.mito.vcf.gz.tbi"

	script:
	"""
	set -euo pipefail

	MITO_CONTIG="${mito_contig}"
	if ! grep -Eq "^>\$MITO_CONTIG([[:space:]]|\$)" "$ref_fa"; then
	    if grep -Eq '^>MT([[:space:]]|\$)' "$ref_fa"; then
	        MITO_CONTIG="MT"
	    elif grep -Eq '^>chrM([[:space:]]|\$)' "$ref_fa"; then
	        MITO_CONTIG="chrM"
	    fi
	fi

	gatk Mutect2 \
	    -R "$ref_fa" \
	    -I "$bam" \
	    -L "\$MITO_CONTIG" \
	    --mitochondria-mode true \
	    -O "${out_prefix}.mito.raw.vcf.gz"

	gatk FilterMutectCalls \
	    -R "$ref_fa" \
	    -V "${out_prefix}.mito.raw.vcf.gz" \
	    -O "${out_prefix}.mito.filtered.vcf.gz"

	bcftools norm -m -any "${out_prefix}.mito.filtered.vcf.gz" -Oz -o "${out_prefix}.mito.vcf.gz"
	tabix -f -p vcf "${out_prefix}.mito.vcf.gz"
	"""
}
