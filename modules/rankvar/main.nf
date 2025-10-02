
process RankVar {
        container ='beoungl/docker_test:rankvar'
	
	input:
	path vcf 
	val phen2gene
	path hpo
	val hpo_directory
	val out_prefix
	val gnomad
	val gq
	val ad
	path input_directory_path
	val input_directory

	output:
	val "${out_prefix}.rank_var.tsv"

	script:
	"""
	source /conda/etc/profile.d/conda.sh
	conda activate rankvar

	python /opt/RankVar/RankVar.py --annovar $input_directory/$vcf --output $input_directory/${out_prefix}_rankvar --hpo_ids $hpo_directory/$hpo --phen2gene $phen2gene --gq $gq --ad $ad --gnomad $gnomad
	
	mv $input_directory/${out_prefix}_rankvar/rank_var.tsv $input_directory/${out_prefix}.rank_var.tsv

	"""
}



