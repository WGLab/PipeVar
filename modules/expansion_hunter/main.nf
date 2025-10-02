
process ExpansionHunter {
	container ='beoungl/docker_test:eh'


        input:
        path bam
        val input_directory
	val out_prefix
	path ref_fa
	val ref_fa_directory
	path output_directory
	val output_directory_full


	script:

	"""
	ExpansionHunter --reads $input_directory/$bam --reference $ref_fa_directory/$ref_fa --variant-catalog /hg38/variant_catalog.json --output-prefix $output_directory_full/$out_prefix
	"""


}



