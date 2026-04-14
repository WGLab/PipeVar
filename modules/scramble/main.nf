// Call mobile-element insertions from short-read BAM alignments using SCRAMBLE.
process SCRAMBLE {
	container = 'beoungl/docker_test:scramble_0.1'

	input:
	tuple path(bam), path(index)
	val out_prefix
	tuple path(ref_fa), path(fa_index)

	output:
	path "${out_prefix}_scramble.vcf", emit: vcf
	path "${out_prefix}_scramble.clusters.txt", optional: true, emit: clusters
	path "${out_prefix}_scramble.MEIs.txt", optional: true, emit: meis

	script:
	def args = task.ext.args ?: ''
	"""
	SCRAMBLE_ROOT=/opt/scramble
	SCRAMBLE_CLUSTER_IDENTIFIER=\$(command -v cluster_identifier || true)
	SCRAMBLE_R=\${SCRAMBLE_ROOT}/cluster_analysis/bin/SCRAMble.R
	SCRAMBLE_INSTALL_DIR=\${SCRAMBLE_ROOT}/cluster_analysis/bin
	SCRAMBLE_MEI_REFS=\${SCRAMBLE_ROOT}/cluster_analysis/resources/MEI_consensus_seqs.fa
	CLUSTERS_OUT=${out_prefix}_scramble.clusters.txt

	if [[ "${bam}" != *.bam ]]; then
	    echo "SCRAMBLE currently supports BAM input only in PipeVar_mito v1" >&2
	    exit 1
	fi
	if [[ -z "\$SCRAMBLE_CLUSTER_IDENTIFIER" ]]; then
	    echo "SCRAMBLE cluster_identifier binary not found in PATH" >&2
	    exit 1
	fi
	if [[ ! -s "\$SCRAMBLE_R" ]]; then
	    echo "SCRAMBLE R script not found at \$SCRAMBLE_R" >&2
	    exit 1
	fi
	if [[ ! -s "\$SCRAMBLE_MEI_REFS" ]]; then
	    echo "SCRAMBLE MEI reference FASTA not found at \$SCRAMBLE_MEI_REFS" >&2
	    exit 1
	fi

	"\$SCRAMBLE_CLUSTER_IDENTIFIER" $bam > "\$CLUSTERS_OUT"

	Rscript --vanilla "\$SCRAMBLE_R" \\
	    --out-name ${out_prefix}_scramble \\
	    --cluster-file "\$CLUSTERS_OUT" \\
	    --install-dir "\$SCRAMBLE_INSTALL_DIR" \\
	    --mei-refs "\$SCRAMBLE_MEI_REFS" \\
	    --ref $ref_fa \\
	    --eval-meis \\
	    $args

	mei_txt=\$(find . -maxdepth 3 -type f -name '*_MEIs.txt' | head -n 1)
	if [[ -n "\$mei_txt" && "\$mei_txt" != "${out_prefix}_scramble.MEIs.txt" ]]; then
	    cp "\$mei_txt" ${out_prefix}_scramble.MEIs.txt
	fi

	vcf_candidate=\$(find . -maxdepth 3 -type f \\( -name "${out_prefix}_scramble*.vcf" -o -name "${out_prefix}_scramble*.vcf.gz" -o -name '*.vcf' -o -name '*.vcf.gz' \\) | head -n 1)
	if [[ -z "\$vcf_candidate" ]]; then
	    echo "SCRAMBLE did not emit a VCF output" >&2
	    exit 1
	fi

	if [[ "\$vcf_candidate" == *.gz ]]; then
	    gunzip -c "\$vcf_candidate" > ${out_prefix}_scramble.vcf
	else
	    cp "\$vcf_candidate" ${out_prefix}_scramble.vcf
	fi
	"""
}
