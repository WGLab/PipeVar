
process ExpansionHunter {
	container ='beoungl/docker_test:eh'


        input:
        tuple path(bam), path(index)
	val out_prefix
	tuple path(ref_fa), path(fa_index)

	output:
	path "${out_prefix}.json"

	script:
	def args = task.ext.args ?: ''

	"""
	ExpansionHunter --reads $bam --reference $ref_fa --variant-catalog $args --output-prefix $out_prefix

		
	
	"""


}



