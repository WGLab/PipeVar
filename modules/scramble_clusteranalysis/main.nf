// Stage 2 SCRAMBLE wrapper: analyze clusters into MEI VCF and text side outputs.
process scramble_clusteranalysis {
	container = 'quay.io/biocontainers/scramble:1.0.2--1'

	input:
	tuple val(meta), path(clusters)
	tuple val(meta2), path(ref_fa), path(fa_index)

	output:
	tuple val(meta), path("${meta.id}_scramble.vcf"), emit: vcf
	tuple val(meta), path("${meta.id}_scramble.MEIs.txt"), optional: true, emit: meis
	tuple val(meta), path("${meta.id}_scramble.PredictedDeletions.txt"), optional: true, emit: dels
	path "versions.yml", emit: versions

	script:
	def args = task.ext.args ?: ''
	def meiRef = task.ext.mei_ref ?: '/opt/scramble/cluster_analysis/resources/MEI_consensus_seqs.fa'
	"""
	SCRAMBLE_R=/opt/scramble/cluster_analysis/bin/SCRAMble.R

	if [[ ! -s "\$SCRAMBLE_R" ]]; then
	    echo "SCRAMBLE R script not found at \$SCRAMBLE_R" >&2
	    exit 1
	fi
	if [[ ! -s "${meiRef}" ]]; then
	    echo "SCRAMBLE MEI reference FASTA not found at ${meiRef}" >&2
	    exit 1
	fi

	Rscript --vanilla "\$SCRAMBLE_R" \\
	    --out-name ${meta.id}_scramble \\
	    --cluster-file $clusters \\
	    --install-dir /opt/scramble/cluster_analysis/bin \\
	    --mei-refs ${meiRef} \\
	    --ref $ref_fa \\
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
