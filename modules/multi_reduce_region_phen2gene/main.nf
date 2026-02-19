
// Batch conversion of top phenotype genes into merged BED intervals.
process multi_phen2gene_filter {
        container ='beoungl/docker_test:phen2gene_filter_0.1.1'
	
	input:
	tuple val(out_prefix), path(phen2gene)
	tuple path(ref_fa), path(fa_index)
	val phen2gene_top_n

	output:
	tuple val(out_prefix),path("${out_prefix}.phen2gene.bed")

	script:
	"""

	bash /gene_to_bed.sh $phen2gene /gtf_file/gencode.v49.annotation.gtf $out_prefix $phen2gene_top_n



	"""
}


