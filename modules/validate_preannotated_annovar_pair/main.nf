// Validate pre-annotated ANNOVAR TXT/VCF pairs before SNP prioritization.
process validate_preannotated_annovar_pair {
	container = 'beoungl/docker_test:validate_preannotated_annovar_pair_0.1'

	input:
	tuple val(out_prefix), path(annovar_txt), path(annovar_vcf)

	output:
	tuple val(out_prefix), path("${out_prefix}.validated.hg38_multianno.txt"), path("${out_prefix}.validated.hg38_multianno.vcf")

	script:
	"""
	validate_preannotated_annovar_pair.py \
	  --sample $out_prefix \
	  --annovar-txt $annovar_txt \
	  --annovar-vcf $annovar_vcf \
	  --validated-txt ${out_prefix}.validated.hg38_multianno.txt \
	  --validated-vcf ${out_prefix}.validated.hg38_multianno.vcf
	"""
}
