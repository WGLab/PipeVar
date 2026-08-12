
// Annotate SV VCF with ANNOVAR and retain exonic records.
process ANNOVAR_SV {
	container ='beoungl/docker_test:truvari_0.5'


	input:
	path vcf
	val out_prefix
	val sv_annotation_mode

	output:
	path "${out_prefix}.exonic.vcf"

	script:
	"""
	source /conda/etc/profile.d/conda.sh
	conda activate truvari



	if [[ "$sv_annotation_mode" == "preannotated" ]]; then
		cp $vcf ${out_prefix}_sv.hg38_multianno.vcf
	else
		perl /annovar/table_annovar.pl $vcf /annovar/humandb/ -buildver hg38 -out ${out_prefix}_sv -remove -protocol refGene -operation gx -nastring . -vcfinput -polish
	fi

	bcftools view -i 'INFO/Func.refGene="exonic"' -Ov -o ${out_prefix}.exonic.vcf ${out_prefix}_sv.hg38_multianno.vcf


	"""

	


}
