
// Run RankVar to prioritize small variants with phenotype-aware evidence.
process RankVar {
        container ='beoungl/docker_test:rankvar'
	
	input:
	path vcf 
	path phen2gene
	path hpo
	val out_prefix
	val gnomad
	val gq
	val ad

	output:
	path "${out_prefix}.rank_var.tsv"

	script:
	"""
	source /conda/etc/profile.d/conda.sh
	conda activate rankvar

	head -n 1 $vcf > ${out_prefix}.rankvar_temp.txt

	grep -Ewi 'exonic|splicing' $vcf >> ${out_prefix}.rankvar_temp.txt

	python /opt/RankVar/RankVar.py --annovar ${out_prefix}.rankvar_temp.txt --output ${out_prefix}_rankvar --hpo_ids $hpo --phen2gene $phen2gene --gq $gq --ad $ad --gnomad $gnomad
	
	#mv ${out_prefix}_rankvar/rank_var.tsv ${out_prefix}.rank_var.tsv
	#This pathogenicity score can change, so is the phen2gene score.
	awk -F'\t' 'NR==1 || \$12 > 0.05' ${out_prefix}_rankvar/rank_var.tsv > ${out_prefix}.rank_var.tsv

	rm ${out_prefix}.rankvar_temp.txt
	"""
}



