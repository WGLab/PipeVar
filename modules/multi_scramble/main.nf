// Combined batch SCRAMBLE wrapper: identify soft-clipped clusters, then analyze them into MEI calls.
process multi_scramble {
	container = 'beoungl/docker_test:scramble_0.1'

	input:
	tuple val(meta), path(bam), path(index)
	tuple val(ref_meta), path(ref_fa), path(fa_index)

	output:
	tuple val(meta), path("${meta.id}_scramble.vcf"), emit: vcf
	tuple val(meta), path("${meta.id}_scramble.clusters.txt"), emit: clusters
	tuple val(meta), path("${meta.id}_scramble.MEIs.txt"), optional: true, emit: meis
	tuple val(meta), path("${meta.id}_scramble.PredictedDeletions.txt"), optional: true, emit: dels
	path "versions.yml", emit: versions

	script:
	def args = task.ext.args ?: ''
	def meiRef = task.ext.mei_ref ?: '/app/cluster_analysis/resources/MEI_consensus_seqs.fa'
	"""
	SCRAMBLE_CLUSTER_IDENTIFIER=\$(command -v cluster_identifier || true)
	SCRAMBLE_R=/app/cluster_analysis/bin/SCRAMble.R

	if [[ -z "\$SCRAMBLE_CLUSTER_IDENTIFIER" ]]; then
	    echo "SCRAMBLE cluster_identifier binary not found in PATH" >&2
	    exit 1
	fi
	if [[ ! -s "\$SCRAMBLE_R" ]]; then
	    echo "SCRAMBLE R script not found at \$SCRAMBLE_R" >&2
	    exit 1
	fi
	if [[ ! -s "${meiRef}" ]]; then
	    echo "SCRAMBLE MEI reference FASTA not found at ${meiRef}" >&2
	    exit 1
	fi

	"\$SCRAMBLE_CLUSTER_IDENTIFIER" "$bam" > ${meta.id}_scramble.clusters.txt

	Rscript --vanilla "\$SCRAMBLE_R" \\
	    --out-name ${meta.id}_scramble \\
	    --cluster-file ${meta.id}_scramble.clusters.txt \\
	    --install-dir /app/cluster_analysis/bin \\
	    --mei-refs ${meiRef} \\
	    --ref $ref_fa \\
	    --eval-dels \\
	    --eval-meis \\
	    $args

	if [[ -f ${meta.id}_scramble_MEIs.txt ]]; then
	    mv ${meta.id}_scramble_MEIs.txt ${meta.id}_scramble.MEIs.txt
	fi
	if [[ -f ${meta.id}_scramble_PredictedDeletions.txt ]]; then
	    mv ${meta.id}_scramble_PredictedDeletions.txt ${meta.id}_scramble.PredictedDeletions.txt
	fi

	cat <<-END_VERSIONS > versions.yml
	"${task.process}":
	  scramble: "unknown"
	END_VERSIONS
	"""
}
