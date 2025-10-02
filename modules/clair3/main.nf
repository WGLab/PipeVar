
process clair3 {
	container ='hkubal/clair3:v1.2.0'


        input:
        path bam
	val input_directory
        val out_prefix
	path ref_fa
	val ref_fa_directory
	path output_directory
	val output_directory_full

	output:
	path "$output_directory/${out_prefix}.clair3.vcf.gz"

	script:
	def args   = task.ext.args ?: ''

	"""
        source /opt/conda/etc/profile.d/conda.sh
	conda activate clair3

	run_clair3.sh --bam_fn=$input_directory/$bam --ref_fn=$ref_fa_directory/$ref_fa --threads=4 $args --output=$output_directory_full/${out_prefix}_clair3
	
	mv $output_directory_full/${out_prefix}_clair3/merge_output.vcf.gz $output_directory_full/${out_prefix}_clair3/${out_prefix}.clair3.vcf.gz

	mv $output_directory_full/${out_prefix}_clair3/${out_prefix}.clair3.vcf.gz $output_directory_full


	"""


}



