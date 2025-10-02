
process sniffles {
        container ='beoungl/docker_test:sniffles'


        input:
        path bam
	val input_directory
        val out_prefix
	path ref_fa
	val ref_fa_directory
	path output_directory
	val output_directory_full
	
	output:
	val "$output_directory_full/${out_prefix}.sniffles.vcf"

	script: 
	"""

	source /conda/etc/profile.d/conda.sh
	conda activate sniffles
	sniffles --allow-overwrite --input $input_directory/$bam --vcf $output_directory_full/${out_prefix}.sniffles.vcf --reference $ref_fa_directory/$ref_fa

	"""

}



