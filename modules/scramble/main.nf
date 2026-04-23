// Combined SCRAMBLE wrapper: identify soft-clipped clusters, then analyze them into MEI calls.
process scramble {
	container = 'beoungl/docker_test:scramble_0.3'

	input:
	tuple val(meta), path(bam), path(index)
	tuple val(ref_meta), path(scramble_ref_dir)

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
	SCRAMBLE_REF_DIR="${scramble_ref_dir}"
	SCRAMBLE_REF="\${SCRAMBLE_REF_DIR}/reference.fa"

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
	if [[ ! -s "\$SCRAMBLE_REF" ]]; then
	    echo "Prepared SCRAMBLE reference FASTA not found at \$SCRAMBLE_REF" >&2
	    exit 1
	fi
	if [[ ! -s "\${SCRAMBLE_REF}.fai" ]]; then
	    echo "Prepared SCRAMBLE FASTA index not found at \${SCRAMBLE_REF}.fai" >&2
	    exit 1
	fi
	if ! command -v blastdbcmd >/dev/null 2>&1; then
	    echo "blastdbcmd not found in PATH; SCRAMBLE deletion validation cannot continue" >&2
	    exit 1
	fi
	if ! blastdbcmd -db "\$SCRAMBLE_REF" -info >/dev/null 2>&1; then
	    echo "Prepared SCRAMBLE BLAST database not found or invalid for prefix \$SCRAMBLE_REF" >&2
	    exit 1
	fi

	if [[ "$bam" == *.cram ]]; then
	    echo "WARNING: SCRAMBLE CRAM support remains provisional and depends on upstream cluster_identifier reference resolution." >&2
	fi

	"\$SCRAMBLE_CLUSTER_IDENTIFIER" "$bam" > ${meta.id}_scramble.clusters.txt

	set +e
	Rscript --vanilla "\$SCRAMBLE_R" \\
	    --out-name ${meta.id}_scramble \\
	    --cluster-file ${meta.id}_scramble.clusters.txt \\
	    --install-dir /app/cluster_analysis/bin \\
	    --mei-refs ${meiRef} \\
	    --ref "\$SCRAMBLE_REF" \\
	    --eval-dels \\
	    --eval-meis \\
	    $args
	scramble_status=\$?
	set -e

	if [[ \$scramble_status -ne 0 && ! -f ${meta.id}_scramble.vcf ]]; then
	    echo "SCRAMBLE R step failed with exit code \$scramble_status" >&2
	    exit \$scramble_status
	fi

	if [[ ! -f ${meta.id}_scramble.vcf ]]; then
	    {
	        echo '##fileformat=VCFv4.2'
	        echo "##reference=\$SCRAMBLE_REF"
	        echo -e '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO'
	    } > ${meta.id}_scramble.vcf
	fi

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
