
// Genotype repeat expansion loci with ExpansionHunter.
process ExpansionHunter {
	container ='beoungl/docker_test:eh'


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


