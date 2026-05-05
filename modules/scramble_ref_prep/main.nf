// Prepare a task-local FASTA plus BLAST database bundle for SCRAMBLE.
process scramble_ref_prep {
	container = 'beoungl/docker_test:scramble_0.3'

	input:
	tuple val(meta), path(ref_fa), path(fa_index)

	output:
	tuple val(meta), path("scramble_ref_db"), emit: ref
	path "versions.yml", emit: versions

	script:
	"""
	REF_DIR="scramble_ref_db"
	REF_BASENAME="reference.fa"
	REF_PATH="\${REF_DIR}/\${REF_BASENAME}"

	mkdir -p "\$REF_DIR"
	ln -s "$ref_fa" "\$REF_PATH"
	ln -s "$fa_index" "\${REF_PATH}.fai"

	makeblastdb -in "\$REF_PATH" -dbtype nucl -parse_seqids -out "\$REF_PATH"
	blastdbcmd -db "\$REF_PATH" -info >/dev/null

	cat <<-END_VERSIONS > versions.yml
	"${task.process}":
	  blast_plus: "\$(makeblastdb -version | head -n 1 | sed 's/^makeblastdb: //')"
	END_VERSIONS
	"""
}
