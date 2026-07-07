
// Annotate SV VCF with ANNOVAR and filter to phenotype-relevant genes.
process ANNOVAR_SV {
	container ='beoungl/docker_test:truvari_0.5'


        input:
        path vcf
	val out_prefix
	path phen2gene
	val bed_file
	val sv_annotation_mode

	output:
	path "${out_prefix}.exonic.vcf"

	script:
	def bed_arg = (bed_file != "null") ? "-bedfile ${bed_file}" : ""
	"""
	source /conda/etc/profile.d/conda.sh
	conda activate truvari



	if [[ "$sv_annotation_mode" == "preannotated" ]]; then
		cp $vcf ${out_prefix}_sv.hg38_multianno.vcf
	else
		perl /annovar/table_annovar.pl $vcf /annovar/humandb/ -buildver hg38 -out ${out_prefix}_sv -remove -protocol refGene -operation gx -nastring . -vcfinput -polish $bed_arg
	fi

	#Filter with Phen2gene score here
	bash /phen2gene_filter.sh $phen2gene ${out_prefix}_sv.hg38_multianno.vcf $out_prefix


	cp ${out_prefix}_sv.phen2gene.vcf ${out_prefix}.exonic.vcf


	"""

	


}
