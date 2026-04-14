
// Batch SV annotation and phenotype-gene filtering with ANNOVAR.
process multi_annovar_sv {
	container ='beoungl/docker_test:truvari_0.1'


        input:
	tuple val(out_prefix), path(vcf), path(phen2gene), val(bed_file)

	output:
	tuple val(out_prefix), path("${out_prefix}.exonic.vcf")


	script:
	def bed_arg = (bed_file != "null") ? "-bedfile ${bed_file}" : ""
	"""
	source /conda/etc/profile.d/conda.sh
	conda activate truvari



        perl /annovar/table_annovar.pl $vcf /annovar/humandb/ -buildver hg38 -out ${out_prefix}_sv -remove -protocol refGene -operation gx -nastring . -vcfinput -polish $bed_arg

	bash /phen2gene_filter.sh $phen2gene ${out_prefix}_sv.hg38_multianno.vcf $out_prefix

	cp ${out_prefix}_sv.phen2gene.vcf ${out_prefix}.exonic.vcf



	"""

	


}
