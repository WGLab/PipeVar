

process SURVIVOR {
        container ='beoungl/docker_test:survivor'

        input:
        path vcf
        val out_prefix
	path input_directory_path
	val input_directory	

	output:
	val "$input_directory/${out_prefix}.bed"


	script:

	"""


	
	SURVIVOR vcftobed $input_directory/$vcf 0 -1 $input_directory/${out_prefix}.int.bed
	awk -F'\t' -v OFS='\t' '{print \$1,\$2,\$5,\$7,\$11}' $input_directory/${out_prefix}.int.bed |
	sed -e 's/\\bINS\\b/insertion/g' \
    -e 's/\\bDEL\\b/deletion/g' \
    -e 's/\\bINV\\b/inversion/g' \
    -e 's/\\bDUP\\b/duplication/g' -e 's/\\bBND\\b/translocation/g' | grep -v 'TRA'  > $input_directory/${out_prefix}.bed

	rm $input_directory/${out_prefix}.int.bed

	"""

}
