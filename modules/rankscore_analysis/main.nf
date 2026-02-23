
// Filter annotated variants by RankScore and emit ClinVar/RankScore subsets.
process Rankscore_analysis {
	container ='beoungl/docker_test:rankscore_0.2.11'


	input:
	path annovar_output
	path phen2_gene
	val out_prefix
	val gnomad_af
	val rankscore_filter
	val gq
	val phen2gene_top_n

	output:
	path "${out_prefix}.rankscore_filtered.tsv" , emit: rankscore
	path "${out_prefix}.clinvar.txt" , emit: clinvar

	"""

	sh /rankscore/clinvar.sh $annovar_output $phen2_gene $out_prefix $gnomad_af $rankscore_filter $phen2gene_top_n $gq

	"""
}
