
process multi_expansionhunter {
	container ='beoungl/docker_test:eh'


        input:
	tuple val(out_prefix), path(bam), path(index_file)
	tuple path(ref_fa), path(fa_index)

	output:
	tuple val(out_prefix), path("${out_prefix}.json")

	script:

	"""
	ExpansionHunter --reads $bam --reference $ref_fa --variant-catalog /hg38/variant_catalog.json --output-prefix $out_prefix

		
	
	"""


}



