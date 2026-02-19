

// Batch conversion/aggregation of SV calls into BED format via SURVIVOR.
process multi_survivor {
        container ='beoungl/docker_test:survivor'

        input:
        tuple val(out_prefix), path(vcf)

	output:
	tuple val(out_prefix), path("${out_prefix}.bed")


	script:

	"""


	
	SURVIVOR vcftobed $vcf 0 -1 ${out_prefix}.int.bed
	awk -F'\t' -v OFS='\t' '{print \$1,\$2,\$5,\$7,\$11}' ${out_prefix}.int.bed |
	sed -e 's/\\bINS\\b/insertion/g' \
    -e 's/\\bDEL\\b/deletion/g' \
    -e 's/\\bINV\\b/inversion/g' \
    -e 's/\\bDUP\\b/duplication/g' -e 's/\\bBND\\b/translocation/g' | grep -v 'TRA' | grep -v 'INS'  > ${out_prefix}.bed

	rm ${out_prefix}.int.bed

	"""

}
