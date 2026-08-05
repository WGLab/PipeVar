

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

	# Always emit a BED output. If no SV rows pass filters, keep header-only file.
	{
	    echo -e "#CHROM\tSTART\tEND\tID\tSVTYPE"
	    awk -F'\t' -v OFS='\t' '{print \$1,\$2,\$5,\$7,\$11}' ${out_prefix}.int.bed | \
	    sed -e 's/\\bINS\\b/insertion/g' \
	        -e 's/\\bDEL\\b/deletion/g' \
	        -e 's/\\bINV\\b/inversion/g' \
	        -e 's/\\bDUP\\b/duplication/g' \
	        -e 's/\\bBND\\b/translocation/g' | \
	    awk 'BEGIN{IGNORECASE=1} \$0 !~ /(^|[[:space:]])TRA([[:space:]]|\$)/ {print}'
	} > ${out_prefix}.bed

	rm -f ${out_prefix}.int.bed

	"""

}
