
// Batch short-read SNP/indel calling with DeepVariant.
process multi_deepvariant {
	container "google/deepvariant:${params.deepvariant_version}${params.GPU?.toString()?.trim()?.toLowerCase() == 'yes' ? '-gpu' : ''}"

        input:
	tuple val(out_prefix), path(bam), path(index_file), path(bed_file)
	tuple path(ref_fa), path(fa_index)

	output:
	tuple val(out_prefix), path("${out_prefix}.deepvariant.vcf.gz")

	script:
	def regions = bed_file ? "--regions=${bed_file}" : ""

	"""

	/opt/deepvariant/bin/run_deepvariant --model_type=WGS --ref=$ref_fa --reads=$bam --output_vcf=${out_prefix}.deepvariant.vcf.gz $regions --num_shards=${task.cpus}



	"""


}
