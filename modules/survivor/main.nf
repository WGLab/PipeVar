

// Convert/normalize SV calls into BED representation with SURVIVOR.
process SURVIVOR {
        container ='beoungl/docker_test:survivor'

        input:
        path vcf
        val out_prefix

	output:
	path "${out_prefix}.bed"


	script:

	"""


	
	SURVIVOR vcftobed $vcf 0 -1 ${out_prefix}.int.bed

	# Always emit a BED output. If no SV rows pass filters, keep header-only file.
	{
	    echo -e "#chrom\tstart\tend\tsvtype\tgene"
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
