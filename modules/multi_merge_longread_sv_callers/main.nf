// Batch merge of long-read raw SV VCFs keyed by sample prefix.
process multi_merge_longread_sv_callers {
	container = 'beoungl/docker_test:longphase_0.2.17'

	input:
	tuple val(out_prefix), path(sv_vcfs)

	output:
	tuple val(out_prefix), path("${out_prefix}.longread_sv.merged.vcf")

	script:
	"""
	vcf_files=( $sv_vcfs )
	if [[ \${#vcf_files[@]} -eq 0 ]]; then
	    echo "No long-read SV VCFs provided for merging" >&2
	    exit 1
	fi

	first_vcf="\${vcf_files[0]}"
	{
		grep '^##' "\$first_vcf"
		for vcf in "\${vcf_files[@]:1}"; do
			grep '^##' "\$vcf" | grep -v '^##fileformat=' || true
		done
		grep '^#CHROM' "\$first_vcf" | tail -n 1
		for vcf in "\${vcf_files[@]}"; do
			grep -v '^#' "\$vcf"
		done
	} | awk 'BEGIN{OFS="\\t"} /^#/ {print; next} NF >= 8 { if (\$3 == ".") \$3 = "SV_" NR; print }' > ${out_prefix}.longread_sv.merged.vcf
	"""
}
