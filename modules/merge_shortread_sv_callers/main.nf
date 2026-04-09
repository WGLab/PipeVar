// Merge Manta and CNVnator SV/CNV calls into a single VCF for downstream annotation.
process merge_shortread_sv_callers {
	container = 'beoungl/docker_test:longphase_0.2.17'

	input:
	path manta_vcf
	path cnvnator_vcf
	val out_prefix

	output:
	path "${out_prefix}.shortread_sv.merged.vcf"

	script:
	"""
	{
		grep '^##' $manta_vcf
		grep '^##' $cnvnator_vcf | grep -v '^##fileformat=' || true
		grep '^#CHROM' $manta_vcf | tail -n 1
		grep -v '^#' $manta_vcf
		grep -v '^#' $cnvnator_vcf
	} | awk 'BEGIN{OFS="\\t"} /^#/ {print; next} NF >= 8 { if (\$3 == ".") \$3 = "SV_" NR; print }' > ${out_prefix}.shortread_sv.merged.vcf
	"""
}
