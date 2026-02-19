
// Convert top phenotype-ranked genes into BED regions for targeted calling.
process phen2gene_filter {
        container ='beoungl/docker_test:phen2gene_filter_0.1.1'
	
	input:
	path phen2gene
	tuple path(ref_fa), path(fa_index)
	val out_prefix
	val phen2gene_top_n

	output:
	path "${out_prefix}.phen2gene.bed"

	script:
	"""

	bash /gene_to_bed.sh $phen2gene /gtf_file/gencode.v49.annotation.gtf $out_prefix $phen2gene_top_n



	"""
}


