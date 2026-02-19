
process multi_sv_prio {
	container='beoungl/docker_test:longphase_0.2.7'

	input:
	tuple val(out_prefix), path(sv_pathogenic), path(annovar_sv_vcf)


	output:
	path "${out_prefix}.prio.vcf"
	
	script:


	"""
	#Final process to organize the results and show them in ACMG guideline format; might take long to develop.

	bash /phenosv_vcf_and_tsv.sh $sv_pathogenic $annovar_sv_vcf $out_prefix


	#Get SNP from RankScore + RankVar + ClinVar, and based on score + homozygosity, rank them.

	python3 /assign_dom_or_rec_sv_only.py ${out_prefix}.phenosv.vcf OMIM ${out_prefix}.assigned.vcf

	python3 /prio_gene_only.py ${out_prefix}.assigned.vcf ${out_prefix}.prio.vcf
	

	#Potential homozygotes are ranked higher than single heterozygotes, but nothing conclusive.
	#Also gene-based.

	#Use similar strategy as the all long-phase analysis as well here.


	"""

}

