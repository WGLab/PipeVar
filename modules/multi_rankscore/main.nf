
// Batch RankScore filtering and ClinVar extraction from ANNOVAR outputs.
process multi_rankscore {
	container ='beoungl/docker_test:rankscore_0.2.15'


	input:
	tuple val(out_prefix), path(annovar_output), path(phen2_gene)
	val gnomad_af
	val rankscore_filter
	val gq
	val phen2gene_top_n

	output:
	tuple val(out_prefix), path("${out_prefix}.rankscore_filtered.tsv"), path("${out_prefix}.clinvar.txt")

	"""


	sh /rankscore/clinvar.sh $annovar_output $phen2_gene $out_prefix $gnomad_af $rankscore_filter $phen2gene_top_n $gq

	"""
}
