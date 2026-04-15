// Stage 1 SCRAMBLE wrapper: identify soft-clipped clusters from short-read BAM/CRAM input.
process scramble_clusteridentifier {
	container = 'quay.io/biocontainers/scramble:1.0.2--1'

	input:
	tuple val(meta), path(bam), path(index)
	tuple val(meta2), path(ref_fa), path(fa_index)

	output:
	tuple val(meta), path("${meta.id}_scramble.clusters.txt"), emit: clusters
	path "versions.yml", emit: versions

	script:
	"""
	SCRAMBLE_CLUSTER_IDENTIFIER=\$(command -v cluster_identifier || true)

	if [[ -z "\$SCRAMBLE_CLUSTER_IDENTIFIER" ]]; then
	    echo "SCRAMBLE cluster_identifier binary not found in PATH" >&2
	    exit 1
	fi

	"\$SCRAMBLE_CLUSTER_IDENTIFIER" "$bam" > ${meta.id}_scramble.clusters.txt

	cat <<-END_VERSIONS > versions.yml
	"${task.process}":
	  scramble: "unknown"
	END_VERSIONS
	"""
}
