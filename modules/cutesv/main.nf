
process cuteSV {
	container ='beoungl/docker_test:cutesv'

        input:
        path bam
	val input_directory
        val out_prefix
	path ref_fa
	val ref_fa_directory
	path output_directory
	val output_directory_full

	output:
	val "$output_directory_full/${out_prefix}.cutesv.vcf.gz"

	script:
	def args   = task.ext.args ?: ''
 
	"""
	source /conda/etc/profile.d/conda.sh
	conda activate cutesv 
	cuteSV $args --sample $out_prefix $input_directory/$bam $ref_fa_directory/$ref_fa $output_directory_full/${out_prefix}.cutesv.vcf .
	"""


}



