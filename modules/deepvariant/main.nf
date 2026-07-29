
// Call SNP/indel variants from short-read alignments using DeepVariant.
process deepvariant {
	container "google/deepvariant:${params.deepvariant_version}${params.GPU?.toString()?.trim()?.toLowerCase() == 'yes' ? '-gpu' : ''}"

        input:
	tuple path(bam), path(bam_index)
        val out_prefix
	tuple path(ref_fa), path(fa_index)
	val bed_file

	output:
	path "${out_prefix}.deepvariant.vcf.gz"

	script:

	"""
	
	if [ $bed_file != "null" ]; then
		/opt/deepvariant/bin/run_deepvariant --model_type=WGS --ref=$ref_fa --reads=$bam --output_vcf=${out_prefix}.deepvariant.vcf.gz --regions=$bed_file --num_shards=${task.cpus}
	else
		/opt/deepvariant/bin/run_deepvariant --model_type=WGS --ref=$ref_fa --reads=$bam --output_vcf=${out_prefix}.deepvariant.vcf.gz --num_shards=${task.cpus}
	fi



	"""


}
