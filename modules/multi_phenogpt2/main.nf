// Batch extraction of HPO terms from clinical notes with GPU-backed PhenoGPT2.
process multi_phenogpt2 {
	container 'beoungl/docker_test:phenogpt2_0.2'

	input:
	tuple val(out_prefix), path(clinical_note)

	output:
	tuple val(out_prefix), path("${out_prefix}_phenogpt2_patient_hpo.txt")

	"""
	export PHENOGPT2_BATCH_SIZE="${params.phenogpt2_batch_size}"
	export PHENOGPT2_CHUNK_BATCH_SIZE="${params.phenogpt2_chunk_batch_size}"
	export PHENOGPT2_WC="${params.phenogpt2_wc}"
	export PHENOGPT2_ATTN_IMPLEMENTATION="${params.phenogpt2_attn_implementation}"
	export PHENOGPT2_NEGATION="${params.phenogpt2_negation}"
	export PHENOGPT2_MODEL_DIR="/opt/phenogpt2/models/phenogpt2"
	export PHENOGPT2_NEGATION_MODEL_DIR="/opt/phenogpt2/models/negation"
	export PHENOGPT2_EMBEDDING_MODEL_DIR="/opt/phenogpt2/models/embedding"
	export PHENOGPT2_CACHE_DIR="${params.phenogpt2_cache_host_path ? '/opt/phenogpt2/cache' : '.phenogpt2_cache'}"
	export PHENOGPT2_MODEL_FINGERPRINT="${params.phenogpt2_model_fingerprint}"
	export PHENOGPT2_NEGATION_MODEL_FINGERPRINT="${params.phenogpt2_negation_model_fingerprint}"
	export PHENOGPT2_EMBEDDING_MODEL_FINGERPRINT="${params.phenogpt2_embedding_model_fingerprint}"
	export PHENOGPT2_COMBINED_MODEL_FINGERPRINT="${params.phenogpt2_combined_model_fingerprint}"

	run_phenogpt2_to_hpo "${clinical_note}" "${out_prefix}"
	"""
}
