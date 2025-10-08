
process Rankscore_analysis {
	container ='beoungl/docker_test:rankscore_0.1'


	input:
	path annovar_output
	val phen2_gene
	val out_prefix
	path input_directory_path
	val input_directory

	output:
	val "${out_prefix}.rankscore_filtered.tsv" , emit: rankscore
	val "${out_prefix}.clinvar.txt" , emit: clinvar

	"""

	sh /rankscore/clinvar.sh $input_directory/$annovar_output $phen2_gene $input_directory/$out_prefix

	"""
}
