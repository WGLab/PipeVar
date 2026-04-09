// Normalize BAM/CRAM input to BAM for callers that require BAM input.
process normalize_shortread_alignment {
	container = 'community.wave.seqera.io/library/htslib_samtools:1.23.1--5b6bb4ede7e612e5'

	input:
	tuple path(bam), path(index)
	val out_prefix
	tuple path(ref_fa), path(fa_index)

	output:
	tuple path("${out_prefix}_cnvnator.bam"), path("${out_prefix}_cnvnator.bam.bai")

	script:
	"""
	if [[ "$bam" == *.cram ]]; then
		samtools view -@ ${task.cpus} -b -T $ref_fa -o ${out_prefix}_cnvnator.unsorted.bam $bam
		samtools sort -@ ${task.cpus} -o ${out_prefix}_cnvnator.bam ${out_prefix}_cnvnator.unsorted.bam
		samtools index -@ ${task.cpus} ${out_prefix}_cnvnator.bam
		rm ${out_prefix}_cnvnator.unsorted.bam
	else
		ln -s $bam ${out_prefix}_cnvnator.bam
		if [[ "$index" == *.bai ]]; then
			ln -s $index ${out_prefix}_cnvnator.bam.bai
		else
			samtools index -@ ${task.cpus} ${out_prefix}_cnvnator.bam
		fi
	fi
	"""
}
