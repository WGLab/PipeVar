
// Annotate SV VCF with ANNOVAR and filter to phenotype-relevant genes.
process ANNOVAR_SV {
	container ='beoungl/docker_test:truvari_0.1'


        input:
        path vcf
	val out_prefix
	path phen2gene

	output:
	path "${out_prefix}.exonic.vcf"

	script:
	"""
	source /conda/etc/profile.d/conda.sh
	conda activate truvari



        perl /annovar/table_annovar.pl $vcf /annovar/humandb/ -buildver hg38 -out ${out_prefix}_sv -remove -protocol refGene -operation gx -nastring . -vcfinput -polish

	#Filter with Phen2gene score here
	bash /phen2gene_filter.sh $phen2gene ${out_prefix}_sv.hg38_multianno.vcf $out_prefix


	bcftools view -h ${out_prefix}_sv.hg38_multianno.vcf > ${out_prefix}.exonic.vcf

	grep -wi 'exonic' ${out_prefix}_sv.phen2gene.vcf >> ${out_prefix}.exonic.vcf


	"""

	


}

