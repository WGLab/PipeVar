
process ANNOVAR_SV {
	container ='beoungl/docker_test:truvari'


        input:
        path vcf
	val out_prefix
	val input_directory
	path output_directory
	val output_directory_full

	output:
	val "$input_directory/${out_prefix}.exonic.vcf"

	script:
	"""
	source /conda/etc/profile.d/conda.sh
	conda activate truvari



        perl /annovar/table_annovar.pl $input_directory/$vcf /annovar/humandb/ -buildver hg38 -out $output_directory/${out_prefix}_sv -remove -protocol refGene -operation gx -nastring . -vcfinput -polish

	bcftools view -h $output_directory/${out_prefix}_sv.hg38_multianno.vcf > $output_directory/${out_prefix}.exonic.vcf

	grep -wi 'exonic' $output_directory/${out_prefix}_sv.hg38_multianno.vcf >> $output_directory/${out_prefix}.exonic.vcf



	"""

	


}


