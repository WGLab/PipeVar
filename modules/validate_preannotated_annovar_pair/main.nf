// Validate pre-annotated ANNOVAR TXT/VCF pairs and pass the staged inputs downstream.
process validate_preannotated_annovar_pair {
	container = 'beoungl/docker_test:validate_preannotated_annovar_pair_0.4'

	input:
	tuple val(out_prefix), path(annovar_txt), path(annovar_vcf)

	output:
	tuple val(out_prefix), path(annovar_txt), path(annovar_vcf)

	script:
	"""
	validate_preannotated_annovar_pair.py \
	  --sample $out_prefix \
	  --annovar-txt $annovar_txt \
	  --annovar-vcf $annovar_vcf
	"""

	stub:
	"""
	true
	"""
}
