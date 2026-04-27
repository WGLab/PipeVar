
// Genotype repeat expansion loci with ExpansionHunter.
process ExpansionHunter {
	conda "${moduleDir}/environment.yml"
	container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/14/14e1d96665f934a98e569fc5a6fa237f98d3753eee2b6f60d0aea8ff9d44f406/data' : 'community.wave.seqera.io/library/expansionhunter:5.0.0--389ada7e191a4fba' }"


        input:
        tuple path(bam), path(index)
	val out_prefix
	tuple path(ref_fa), path(fa_index)
	path variant_catalog

	output:
	path "${out_prefix}.json.gz", emit: json
	path "${out_prefix}.vcf.gz", emit: vcf
	path "${out_prefix}_realigned.bam", emit: bam, optional: true

	script:
	def args = task.ext.args ?: ''
	def args2 = task.ext.args2 ?: ''

	"""
	ExpansionHunter \\
	  ${args} \\
	  --reads ${bam} \\
	  --output-prefix ${out_prefix} \\
	  --reference ${ref_fa} \\
	  --variant-catalog ${variant_catalog}
	bgzip --threads ${task.cpus} ${args2} ${out_prefix}.vcf
	bgzip --threads ${task.cpus} ${args2} ${out_prefix}.json
	"""


}


