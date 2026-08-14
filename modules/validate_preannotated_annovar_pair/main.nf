// Validate pre-annotated ANNOVAR TXT/VCF pairs without materializing copies.
process validate_preannotated_annovar_pair_check {
	container = 'beoungl/docker_test:validate_preannotated_annovar_pair_0.4'

	input:
	tuple val(out_prefix), path(annovar_txt), path(annovar_vcf)

	output:
	tuple val(out_prefix), val(true)

	script:
	"""
	validate_preannotated_annovar_pair.py \
	  --sample $out_prefix \
	  --annovar-txt $annovar_txt \
	  --annovar-vcf $annovar_vcf
	"""
}

// Preserve the existing downstream tuple contract while using validation only as a gate.
workflow validate_preannotated_annovar_pair {
	take:
	input_pairs

	main:
	validation_ok = validate_preannotated_annovar_pair_check(input_pairs)
	validated_pairs = input_pairs
		.join(validation_ok, failOnMismatch: true, failOnDuplicate: true)
		.map { out_prefix, annovar_txt, annovar_vcf, ok ->
			tuple(out_prefix, annovar_txt, annovar_vcf)
		}

	emit:
	validated_pairs
}
