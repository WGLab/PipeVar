
// Filter annotated variants by RankScore and emit ClinVar/RankScore subsets.
process Rankscore_analysis {
	container ='beoungl/docker_test:rankscore_0.3.0'


	input:
	path annovar_output
	path phen2_gene
	val out_prefix
	val gnomad_af
	val rankscore_filter
	val rankscore_softwares
	val gq
	val ad
	val phen2gene_top_n

	output:
	path "${out_prefix}.rankscore_filtered.tsv" , emit: rankscore
	path "${out_prefix}.clinvar.txt" , emit: clinvar

	script:
	"""

	bash /rankscore/clinvar.sh $annovar_output $phen2_gene $out_prefix $gnomad_af $rankscore_filter $phen2gene_top_n $gq $ad $rankscore_softwares

	"""
}

// Validate imported ANNOVAR TXT/VCF inputs in the same task that runs RankScore.
// Keeping validation here avoids copying large VCFs through a standalone gate.
process Rankscore_analysis_preannotated {
	container ='beoungl/docker_test:rankscore_0.3.1'

	input:
	path annovar_output
	path annovar_vcf
	path phen2_gene
	val out_prefix
	val gnomad_af
	val rankscore_filter
	val rankscore_softwares
	val gq
	val ad
	val phen2gene_top_n

	output:
	path "${out_prefix}.rankscore_filtered.tsv" , emit: rankscore
	path "${out_prefix}.clinvar.txt" , emit: clinvar

	script:
	"""

	python3 /rankscore/validate_preannotated_annovar_pair.py --sample $out_prefix --annovar-txt $annovar_output --annovar-vcf $annovar_vcf
	bash /rankscore/clinvar.sh $annovar_output $phen2_gene $out_prefix $gnomad_af $rankscore_filter $phen2gene_top_n $gq $ad $rankscore_softwares

	"""

	stub:
	"""
	touch ${out_prefix}.rankscore_filtered.tsv ${out_prefix}.clinvar.txt
	"""
}
