
process truvari {
        container ='beoungl/docker_test:truvari'


        input:
        path cutesv_vcf
        path sniffles_vcf
	path ref_fa
	val ref_fa_directory
	val out_prefix
	path output_directory_full
	val input_directory

	output:
	val "$input_directory/${out_prefix}_truvari.exonic.vcf"

	script:
	"""
	source /conda/etc/profile.d/conda.sh
	conda activate truvari

	bgzip $input_directory/$cutesv_vcf
	bgzip $input_directory/$sniffles_vcf

	tabix $input_directory/${cutesv_vcf}.gz
	tabix $input_directory/${sniffles_vcf}.gz


	truvari bench -b $input_directory/${cutesv_vcf}.gz -c $input_directory/${sniffles_vcf}.gz -o $input_directory/${out_prefix}_truvari --reference $ref_fa_directory/$ref_fa --pctseq 0 --pctsize 0.5

	perl /annovar/table_annovar.pl $input_directory/${out_prefix}_truvari/tp-base.vcf.gz /annovar/humandb/ -buildver hg38 -out $input_directory/${out_prefix}_truvari -remove -protocol refGene -operation gx -nastring . -vcfinput -polish


	bcftools view -h $input_directory/${out_prefix}_truvari.hg38_multianno.vcf > $input_directory/${out_prefix}_truvari.exonic.vcf

	grep -wi 'exonic' $input_directory/${out_prefix}_truvari.hg38_multianno.vcf >> $input_directory/${out_prefix}_truvari.exonic.vcf

	"""


}



