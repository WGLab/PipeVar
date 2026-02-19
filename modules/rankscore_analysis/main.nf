
process Rankscore_analysis {
	container ='beoungl/docker_test:rankscore_0.2.11'


	input:
	path annovar_output
	path phen2_gene
	val out_prefix
	val rankscore_filter
	val gnomad_af

	output:
	path "${out_prefix}.rankscore_filtered.tsv" , emit: rankscore
	path "${out_prefix}.clinvar.txt" , emit: clinvar

	"""

	sh /rankscore/clinvar.sh $annovar_output $phen2_gene $out_prefix $gnomad_af $rankscore_filter

	"""
}
