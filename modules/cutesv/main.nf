
process cuteSV {
	container ='beoungl/docker_test:cutesv'

        input:
	tuple path(bam), path(index)
        val out_prefix
	tuple path(ref_fa), path(fa_index)

	output:
	path "${out_prefix}.cutesv.vcf"

	script:
	def args   = task.ext.args ?: ''
 
	"""
	source /conda/etc/profile.d/conda.sh
	conda activate cutesv 
	cuteSV $args --sample $out_prefix $bam $ref_fa ${out_prefix}.cutesv.vcf .
	"""


}



