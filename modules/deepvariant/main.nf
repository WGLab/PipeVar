
process deepvariant {
	container='google/deepvariant:1.9.0'

        input:
	tuple path(bam), path(index)
        val out_prefix
	tuple path(ref_fa), path(index)
	val bed_file

	output:
	path "${out_prefix}.deepvariant.vcf.gz"

	script:

	"""
	
	if [ $bed_file != "null" ]; then
		/opt/deepvariant/bin/run_deepvariant --model_type=WGS --ref=$ref_fa --reads=$bam --output_vcf=${out_prefix}.deepvariant.vcf.gz --regions=$bed_file
	else
		/opt/deepvariant/bin/run_deepvariant --model_type=WGS --ref=$ref_fa --reads=$bam --output_vcf=${out_prefix}.deepvariant.vcf.gz
	fi



	"""


}



