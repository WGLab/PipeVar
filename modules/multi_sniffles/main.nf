
// Batch structural-variant calling from long-read BAM inputs with Sniffles.
process multi_sniffles {
        container ='beoungl/docker_test:sniffles_0.1'


        input:
	tuple val(out_prefix), path(bam), path(index_file)
	tuple path(ref_fa), path(fa_index)
	
	output:
	tuple val(out_prefix), path("${out_prefix}.sniffles.vcf")

	script: 
	"""

	sniffles --allow-overwrite --output-rnames --input $bam --vcf ${out_prefix}.sniffles.vcf --reference $ref_fa

	"""

}


