
// Compare/merge long-read SV calls with Truvari before downstream SV prioritization.
process truvari {
        container ='beoungl/docker_test:truvari'


        input:
        path cutesv_vcf
        path sniffles_vcf
	tuple path(ref_fa), path(fa_index)
	val out_prefix

	output:
	path "${out_prefix}_truvari.exonic.vcf"

	script:
	"""
	source /conda/etc/profile.d/conda.sh
	conda activate truvari

	bgzip $cutesv_vcf
	bgzip $sniffles_vcf

	tabix ${cutesv_vcf}.gz
	tabix ${sniffles_vcf}.gz


	truvari bench -b ${cutesv_vcf}.gz -c ${sniffles_vcf}.gz -o ${out_prefix}_truvari --reference $ref_fa --pctseq 0 --pctsize 0.5

	perl /annovar/table_annovar.pl ${out_prefix}_truvari/tp-base.vcf.gz /annovar/humandb/ -buildver hg38 -out ${out_prefix}_truvari -remove -protocol refGene -operation gx -nastring . -vcfinput -polish


	bcftools view -h ${out_prefix}_truvari.hg38_multianno.vcf > ${out_prefix}_truvari.exonic.vcf

	grep -wi 'exonic' ${out_prefix}_truvari.hg38_multianno.vcf >> ${out_prefix}_truvari.exonic.vcf

	"""


}



