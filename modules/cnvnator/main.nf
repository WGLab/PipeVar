// Call short-read copy-number variants with CNVnator and emit VCF for PipeVar SV prioritization.
process CNVnator {
	container = 'community.wave.seqera.io/library/cnvnator:0.4.1--5a467cfadbbc668d'

	input:
	tuple path(bam), path(index)
	val out_prefix
	tuple path(ref_fa), path(fa_index)
	val bin_size

	output:
	path "${out_prefix}_cnvnator.root", emit: root
	path "${out_prefix}_cnvnator.tab", emit: tab
	path "${out_prefix}_cnvnator.vcf", emit: vcf

	script:
	def args = task.ext.args ?: ''
	"""
	cnvnator -root ${out_prefix}_cnvnator.root -tree $bam $args
	cnvnator -root ${out_prefix}_cnvnator.root -his $bin_size -fasta $ref_fa $args
	cnvnator -root ${out_prefix}_cnvnator.root -stat $bin_size $args
	cnvnator -root ${out_prefix}_cnvnator.root -partition $bin_size $args
	cnvnator -root ${out_prefix}_cnvnator.root -call $bin_size $args > ${out_prefix}_cnvnator.tab
	cnvnator2VCF.pl ${out_prefix}_cnvnator.tab > ${out_prefix}_cnvnator.vcf
	"""
}
