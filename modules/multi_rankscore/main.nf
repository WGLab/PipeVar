
// Batch RankScore filtering and ClinVar extraction from ANNOVAR outputs.
process multi_rankscore {
	container ='beoungl/docker_test:rankscore_0.3.0'


	input:
	tuple val(out_prefix), path(annovar_output), path(phen2_gene)
	val gnomad_af
	val rankscore_filter
	val rankscore_softwares
	val gq
	val ad
	val phen2gene_top_n

	output:
	tuple val(out_prefix), path("${out_prefix}.rankscore_filtered.tsv"), path("${out_prefix}.clinvar.txt")

	script:
	"""


	bash /rankscore/clinvar.sh $annovar_output $phen2_gene $out_prefix $gnomad_af $rankscore_filter $phen2gene_top_n $gq $ad $rankscore_softwares

	"""
}

// Batch RankScore for imported ANNOVAR pairs. Each sample validates and ranks
// independently, so one slow sample cannot hold back the rest of the cohort.
process multi_rankscore_preannotated {
	container ='beoungl/docker_test:rankscore_0.3.1'

	input:
	tuple val(out_prefix), path(annovar_output), path(annovar_vcf), path(phen2_gene)
	val gnomad_af
	val rankscore_filter
	val rankscore_softwares
	val gq
	val ad
	val phen2gene_top_n

	output:
	tuple val(out_prefix), path("${out_prefix}.rankscore_filtered.tsv"), path("${out_prefix}.clinvar.txt")

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
