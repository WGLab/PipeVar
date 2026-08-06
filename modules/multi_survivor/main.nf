

// Batch conversion/aggregation of SV calls into BED format via SURVIVOR.
process multi_survivor {
        container ='beoungl/docker_test:survivor_0.1'

        input:
        tuple val(out_prefix), path(vcf)

	output:
	tuple val(out_prefix), path("${out_prefix}.canonical.bed"), path("${out_prefix}.phenosv.bed"), path("${out_prefix}.phenosv.bedpe"), path("${out_prefix}.phenosv.members.tsv")


	script:

	"""


	
	SURVIVOR vcftobed $vcf 0 -1 ${out_prefix}.int.bed

	python3 /prepare_phenosv_inputs.py $vcf ${out_prefix}.int.bed $out_prefix

	rm -f ${out_prefix}.int.bed

	"""

}
