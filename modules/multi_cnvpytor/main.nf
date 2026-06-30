// Batch long-read CNV calling with CNVpytor keyed by sample prefix.
process multi_cnvpytor {
	container = 'quay.io/biocontainers/cnvpytor:1.3.2--pyhdfd78af_0'

	input:
	tuple val(out_prefix), path(bam), path(index), path(snp_vcf)
	tuple path(ref_fa), path(fa_index)
	val cnvpytor_reference_genome
	path cnvpytor_reference_conf
	val use_baf
	val bin_sizes
	val primary_bin
	val min_size

	output:
	tuple val(out_prefix), path("${out_prefix}.pytor"), emit: pytor
	tuple val(out_prefix), path("${out_prefix}.cnvpytor.tsv"), emit: tsv
	tuple val(out_prefix), path("${out_prefix}.cnvpytor.vcf"), emit: vcf

	script:
	def args = task.ext.args ?: ''
	def ref_arg = bam.name.toLowerCase().endsWith('.cram') ? "-T ${ref_fa}" : ''
	def rg_arg = cnvpytor_reference_genome ? "-rg ${cnvpytor_reference_genome}" : ''
	def conf_arg = cnvpytor_reference_conf ? "-conf ${cnvpytor_reference_conf}" : ''
	"""
	cnvpytor -root ${out_prefix}.pytor -rd $bam ${ref_arg} ${conf_arg} $args
	if [ -n "${rg_arg}" ]; then
		cnvpytor -root ${out_prefix}.pytor ${rg_arg} ${conf_arg} $args
	fi
	cnvpytor -root ${out_prefix}.pytor -his ${bin_sizes} ${conf_arg} $args
	if [ "$use_baf" = "yes" ] && [ -s "$snp_vcf" ]; then
		vcf_sample=\$( (gzip -dc "$snp_vcf" 2>/dev/null || cat "$snp_vcf") | awk '/^#CHROM/ {print \$10; exit}' )
		if [ -n "\$vcf_sample" ]; then
			cnvpytor -root ${out_prefix}.pytor -snp $snp_vcf -sample "\$vcf_sample" ${conf_arg} $args
			cnvpytor -root ${out_prefix}.pytor -pileup $bam ${ref_arg} ${conf_arg} $args
			cnvpytor -root ${out_prefix}.pytor -baf ${bin_sizes} ${conf_arg} $args
		else
			echo "WARNING: CNVpytor BAF mode requested but no VCF sample column was found for $snp_vcf; continuing with RD-only output." >&2
		fi
	fi
	cnvpytor -root ${out_prefix}.pytor -partition ${bin_sizes} ${conf_arg} $args
	cnvpytor -root ${out_prefix}.pytor -call ${primary_bin} ${conf_arg} $args > ${out_prefix}.cnvpytor.raw.tsv

	cat <<EOF | cnvpytor -root ${out_prefix}.pytor -view ${primary_bin} ${conf_arg} $args
	set Q0_range 0 0.5
	set size_range ${min_size} inf
	set print_filename ${out_prefix}.cnvpytor.tsv
	print calls
	set print_filename ${out_prefix}.cnvpytor.vcf
	print calls
	EOF
	"""
}
