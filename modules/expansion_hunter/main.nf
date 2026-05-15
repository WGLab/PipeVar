
// Genotype repeat expansion loci with ExpansionHunter.
process ExpansionHunter {
	container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/14/14e1d96665f934a98e569fc5a6fa237f98d3753eee2b6f60d0aea8ff9d44f406/data' : 'community.wave.seqera.io/library/expansionhunter:5.0.0--389ada7e191a4fba' }"

        input:
        tuple path(bam), path(index)
	val out_prefix
	tuple path(ref_fa), path(fa_index)
	path variant_catalog

	output:
	path "${out_prefix}.json"

	script:

	"""
	ExpansionHunter --reads "$bam" --reference "$ref_fa" --variant-catalog "$variant_catalog" --output-prefix "$out_prefix"
	"""


}


