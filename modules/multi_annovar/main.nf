
// Batch ANNOVAR annotation for per-sample small-variant VCFs.
process multi_annovar {
	container ='beoungl/docker_test:annovar'

        input:
	tuple val(out_prefix), path(vcf)


	output:
	tuple val(out_prefix), path("${out_prefix}.hg38_multianno.txt"), path("${out_prefix}.hg38_multianno.vcf")

	script:
	"""

        perl /annovar/table_annovar.pl $vcf /annovar/humandb/ -buildver hg38 -out $out_prefix -remove -protocol refGene,cytoBand,exac03,avsnp147,dbnsfp47a,gnomad41_exome,gnomad41_genome,clinvar_20240917,GTEx_v8_eQTL,GTEx_v8_sQTL -operation gx,r,f,f,f,f,f,f,f,f -nastring . -vcfinput -polish -otherinfo




	"""

	


}


