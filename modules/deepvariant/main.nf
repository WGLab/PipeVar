
process deepvariant {
	container='google/deepvariant:1.9.0'

        input:
        path bam
	val input_directory
        val out_prefix
	path ref_fa
	val ref_fa_directory
	path output_directory
	val output_directory_full

	output:
	path "$output_directory/${out_prefix}.deepvariant.vcf.gz"

	script:

	"""

	/opt/deepvariant/bin/run_deepvariant --model_type=WGS --ref=$ref_fa_directory/$ref_fa --reads=$input_directory/$bam --output_vcf=$output_directory_full/${out_prefix}.deepvariant.vcf.gz



	"""


}



