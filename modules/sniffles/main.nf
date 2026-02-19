
// Call structural variants from long-read alignments using Sniffles.
process sniffles {
        container ='beoungl/docker_test:sniffles'


        input:
	tuple path(bam), path(index)
        val out_prefix
	tuple path(ref_fa), path(fa_index)
	
	output:
	path "${out_prefix}.sniffles.vcf"

	script: 
	"""

	source /conda/etc/profile.d/conda.sh
	conda activate sniffles
	sniffles --allow-overwrite --input $bam --vcf ${out_prefix}.sniffles.vcf --reference $ref_fa

	"""

}



