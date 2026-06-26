// Experimental long-read CNV calling with CNVpytor.
// Supports RD-only mode and optional SNP/BAF-aware mode for all-longphase flows.
process CNVpytor {
	container = 'quay.io/biocontainers/cnvpytor:1.3.2--pyhdfd78af_0'

	input:
	tuple path(bam), path(index)
	val out_prefix
	tuple path(ref_fa), path(fa_index)
	path snp_vcf
	val use_baf
	val bin_sizes
	val primary_bin
	val min_size

	output:
	path "${out_prefix}.pytor", emit: pytor
	path "${out_prefix}.cnvpytor.tsv", emit: tsv
	path "${out_prefix}.cnvpytor.vcf", emit: vcf

	script:
	def args = task.ext.args ?: ''
	def ref_arg = bam.name.endsWith('.cram') ? "-T ${ref_fa}" : ''
	"""
	cnvpytor -root ${out_prefix}.pytor -rd $bam ${ref_arg} $args
	cnvpytor -root ${out_prefix}.pytor -his ${bin_sizes} $args
	if [ "$use_baf" = "yes" ] && [ -s "$snp_vcf" ]; then
		vcf_sample=\$( (gzip -dc "$snp_vcf" 2>/dev/null || cat "$snp_vcf") | awk '/^#CHROM/ {print \$10; exit}' )
		if [ -n "\$vcf_sample" ]; then
			cnvpytor -root ${out_prefix}.pytor -snp $snp_vcf -sample "\$vcf_sample" $args
			cnvpytor -root ${out_prefix}.pytor -pileup $bam ${ref_arg} $args
			cnvpytor -root ${out_prefix}.pytor -baf ${bin_sizes} $args
		else
			echo "WARNING: CNVpytor BAF mode requested but no VCF sample column was found for $snp_vcf; continuing with RD-only output." >&2
		fi
	fi
	cnvpytor -root ${out_prefix}.pytor -partition ${bin_sizes} $args
	cnvpytor -root ${out_prefix}.pytor -call ${primary_bin} > ${out_prefix}.cnvpytor.raw.tsv

	cat <<EOF | cnvpytor -root ${out_prefix}.pytor -view ${primary_bin}
	set Q0_range 0 0.5
	set size_range ${min_size} inf
	set print_filename ${out_prefix}.cnvpytor.tsv
	print calls
	set print_filename ${out_prefix}.cnvpytor.vcf
	print calls
	EOF
	"""
}
