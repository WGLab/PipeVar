
process nanocaller {
        container ='beoungl/docker_test:nanocaller'


        input:
        path bam
	val input_directory
        val out_prefix
	path ref_fa
	val ref_fa_directory
	path output_directory
	val output_directory_full


	output:
	path "$output_directory/${out_prefix}.nanocaller.vcf.gz"
	
	script:
	def args   = task.ext.args ?: ''


	"""
	source /conda/etc/profile.d/conda.sh
        conda activate nanocaller

	NanoCaller --bam $input_directory/$bam --ref $ref_fa_directory/$ref_fa --cpu 4 --output $output_directory_full/$out_prefix $args

	mv $output_directory_full/$out_prefix/variant_calls.vcf.gz $output_directory_full/$out_prefix/${out_prefix}.nanocaller.vcf.gz

	mv $output_directory_full/$out_prefix/${out_prefix}.nanocaller.vcf.gz $output_directory_full


	"""


}



