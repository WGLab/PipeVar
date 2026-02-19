
process multi_deepvariant {
	container='google/deepvariant:1.9.0'

        input:
	tuple val(out_prefix), path(bam), path(index_file)
	tuple path(ref_fa), path(fa_index)
	val bed_file

	output:
	tuple val(out_prefix), path("${out_prefix}.deepvariant.vcf.gz")

	script:

	"""

	if [ $bed_file != "null" ]; then
		/opt/deepvariant/bin/run_deepvariant --model_type=WGS --ref=$ref_fa --reads=$bam --output_vcf=${out_prefix}.deepvariant.vcf.gz --regions=$bed_file
	else
		/opt/deepvariant/bin/run_deepvariant --model_type=WGS --ref=$ref_fa --reads=$bam --output_vcf=${out_prefix}.deepvariant.vcf.gz
	fi



	"""


}



