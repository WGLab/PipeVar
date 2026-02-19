
// Batch structural-variant calling from long-read BAM inputs with Sniffles.
process multi_sniffles {
        container ='beoungl/docker_test:sniffles'


        input:
	tuple val(out_prefix), path(bam), path(index_file)
	tuple path(ref_fa), path(fa_index)
	
	output:
	tuple val(out_prefix), path("${out_prefix}.sniffles.vcf")

	script: 
	"""

	source /conda/etc/profile.d/conda.sh
	conda activate sniffles
	sniffles --allow-overwrite --input $bam --vcf ${out_prefix}.sniffles.vcf --reference $ref_fa

	"""

}



