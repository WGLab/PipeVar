// Extract HPO terms from clinical notes with GPU-backed PhenoGPT2.
process phenogpt2 {
	container 'beoungl/docker_test:phenogpt2_0.1'

	input:
	path clinical_note
	val out_prefix

	output:
	path "${out_prefix}_phenogpt2_patient_hpo.txt"

	"""
	export PHENOGPT2_BATCH_SIZE="${params.phenogpt2_batch_size}"
	export PHENOGPT2_CHUNK_BATCH_SIZE="${params.phenogpt2_chunk_batch_size}"
	export PHENOGPT2_WC="${params.phenogpt2_wc}"
	export PHENOGPT2_ATTN_IMPLEMENTATION="${params.phenogpt2_attn_implementation}"
	export PHENOGPT2_NEGATION="${params.phenogpt2_negation}"

	run_phenogpt2_to_hpo $clinical_note ${out_prefix}
	"""
}
