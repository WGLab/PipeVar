
// Batch SNP evidence merging (ClinVar/RankScore/RankVar) to prioritized VCF.
process multi_snp_prio {
	container ='beoungl/docker_test:longphase_0.2.9'

	input:
	tuple val(out_prefix), path(snv_rankscore), path(snv_pathogenic), path(snv_rankvar), path(annovar_vcf), path(hpo_path)


	output:
	path "${out_prefix}.prio.vcf"
	
	script:


	"""
	#Final process to organize the results and show them in ACMG guideline format; might take long to develop.

	bash /clinvar_vcf_and_txt.sh $snv_pathogenic $annovar_vcf $out_prefix

        bash /rankscore_vcf_and_txt.sh $snv_rankscore $annovar_vcf $out_prefix

        bash /rankvar_vcf_and_tsv.sh $snv_rankvar $annovar_vcf $out_prefix


	#Get SNP from RankScore + RankVar + ClinVar, and based on score + homozygosity, rank them.

	python3 /assign_dom_or_rec_snp_only.py ${out_prefix}.clinvar.vcf ${out_prefix}.rankscore.vcf ${out_prefix}.rankvar.vcf $hpo_path OMIM ${out_prefix}.assigned.vcf

	python3 /prio_gene_only.py ${out_prefix}.assigned.vcf ${out_prefix}.prio.vcf
	

	#Potential homozygotes are ranked higher than single heterozygotes, but nothing conclusive.
	#Also gene-based.

	#Use similar strategy as the all long-phase analysis as well here.


	"""

}
