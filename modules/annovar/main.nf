
process ANNOVAR {
	container ='beoungl/docker_test:annovar'

        input:
        path vcf
	val input_directory
	val out_prefix
	path output_directory
	val output_directory_full


	output:
	val "$input_directory/${out_prefix}.hg38_multianno.txt"

	script:
	"""

        perl /annovar/table_annovar.pl $input_directory/$vcf /annovar/humandb/ -buildver hg38 -out $output_directory_full/$out_prefix -remove -protocol refGene,cytoBand,exac03,avsnp147,dbnsfp47a,gnomad41_exome,gnomad41_genome,clinvar_20240917,GTEx_v8_eQTL,GTEx_v8_sQTL -operation gx,r,f,f,f,f,f,f,f,f -nastring . -vcfinput -polish -otherinfo




	"""

	


}


