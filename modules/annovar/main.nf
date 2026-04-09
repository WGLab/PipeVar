
// Annotate small variants with ANNOVAR databases and emit annotated TXT/VCF.
process ANNOVAR {
	container ='beoungl/docker_test:annovar'

        input:
        path vcf
	val out_prefix
	val bed_file


	output:
	path "${out_prefix}.hg38_multianno.txt", emit: txt_output
	path "${out_prefix}.hg38_multianno.vcf", emit: vcf_output

	script:
	def bed_arg = (bed_file != "null") ? "-bedfile ${bed_file}" : ""
	"""

        perl /annovar/table_annovar.pl $vcf /annovar/humandb/ -buildver hg38 -out $out_prefix -remove -protocol refGene,cytoBand,exac03,avsnp147,dbnsfp47a,gnomad41_exome,gnomad41_genome,clinvar_20240917,GTEx_v8_eQTL,GTEx_v8_sQTL -operation gx,r,f,f,f,f,f,f,f,f -nastring . -vcfinput -polish -otherinfo $bed_arg




	"""

	


}

