
process clinvar_snp {
	executor 'slurm'
	cpus 8
	memory '128GB'
	time '24h'
	
	input:
	val vcf 
	val phen2gene
	val out_prefix
	val input_directory


	script:
	"""

	#filter for top 500 phen2gene results, and see if there are any clinvar pathogenic results that are available here.	
	
	head -n 500 $phen2gene | grep -v 'Gene' > $phene2_gene_top_500
	while IFS= read -r line
	do
		#Change this so clinvar column gets grepped for pathogenic only
        	grep -wi "$line" $vcf | grep -v 'Pathogenic' >> ${out_prefix}.clinvar.txt
	done < "$phen2_gene_top_500"

	rm $phen2_gene_top_500	 

	"""
}



