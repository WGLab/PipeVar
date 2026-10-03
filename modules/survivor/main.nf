// Convert/normalize SV calls into BED representation with SURVIVOR.
process survivor {
    container 'beoungl/docker_test:survivor_0.2'

    input:
    tuple val(out_prefix), path(vcf)

    output:
    tuple val(out_prefix), path("${out_prefix}.canonical.bed"), path("${out_prefix}.phenosv.bed"), path("${out_prefix}.phenosv.bedpe"), path("${out_prefix}.phenosv.members.tsv"), emit: phenosv_inputs

    script:

	"""


	
	SURVIVOR vcftobed $vcf 0 -1 ${out_prefix}.int.bed

	python3 /prepare_phenosv_inputs.py $vcf ${out_prefix}.int.bed $out_prefix

	rm -f ${out_prefix}.int.bed

	"""

}
