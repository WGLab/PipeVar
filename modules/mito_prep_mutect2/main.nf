process mito_prep_mutect2 {
	container = 'beoungl/docker_test:mito_mutect2_0.1'

	input:
	tuple path(bam), path(index_file)
	val out_prefix
	tuple path(ref_fa), path(fa_index), path(dict_index), path(bwa_amb), path(bwa_ann), path(bwa_bwt), path(bwa_pac), path(bwa_sa)
	val mito_contig

	output:
	tuple path("${out_prefix}.mito.prepped.bam"), path("${out_prefix}.mito.prepped.bam.bai"), path("${out_prefix}.mito.dup_metrics.txt")

	script:
	"""
	MITO_CONTIG="${mito_contig}"
	INPUT_BAM="$bam"
	INPUT_INDEX="$index_file"

	if [[ "$bam" == *.cram ]]; then
	    # DRAGEN CRAMs must be decoded against the exact original reference bundle.
	    samtools view -@ ${task.cpus} -b -T "$ref_fa" -o "${out_prefix}.mito.input.unsorted.bam" "$bam"
	    samtools sort -@ ${task.cpus} -o "${out_prefix}.mito.input.bam" "${out_prefix}.mito.input.unsorted.bam"
	    samtools index -@ ${task.cpus} "${out_prefix}.mito.input.bam"
	    INPUT_BAM="${out_prefix}.mito.input.bam"
	    INPUT_INDEX="${out_prefix}.mito.input.bam.bai"
	fi

	if ! samtools view -H "\$INPUT_BAM" | grep -q '^@RG'; then
	    gatk AddOrReplaceReadGroups \
	        -I "\$INPUT_BAM" \
	        -O "${out_prefix}.mito.rg.bam" \
	        -RGID "${out_prefix}" \
	        -RGLB "${out_prefix}" \
	        -RGPL "ILLUMINA" \
	        -RGPU "${out_prefix}.unit1" \
	        -RGSM "${out_prefix}" \
	        --CREATE_INDEX true
	    INPUT_BAM="${out_prefix}.mito.rg.bam"
	    INPUT_INDEX="${out_prefix}.mito.rg.bam.bai"
	fi

	if ! samtools idxstats "\$INPUT_BAM" | cut -f1 | grep -qx "\$MITO_CONTIG"; then
	    if samtools idxstats "\$INPUT_BAM" | cut -f1 | grep -qx "MT"; then
	        MITO_CONTIG="MT"
	    elif samtools idxstats "\$INPUT_BAM" | cut -f1 | grep -qx "chrM"; then
	        MITO_CONTIG="chrM"
	    elif samtools idxstats "\$INPUT_BAM" | cut -f1 | grep -qx "M"; then
	        MITO_CONTIG="M"
	    elif samtools idxstats "\$INPUT_BAM" | cut -f1 | grep -qx "chrMT"; then
	        MITO_CONTIG="chrMT"
	    else
	        echo "Unable to find mitochondrial contig in \$INPUT_BAM" >&2
	        exit 1
	    fi
	fi

	gatk PrintReads -R "$ref_fa" -I "\$INPUT_BAM" -L "\$MITO_CONTIG" -O "${out_prefix}.mito.subset.bam"
	samtools index -@ ${task.cpus} "${out_prefix}.mito.subset.bam"
	gatk RevertSam \
	    -I "${out_prefix}.mito.subset.bam" \
	    -O "${out_prefix}.mito.ubam" \
	    --REMOVE_ALIGNMENT_INFORMATION true \
	    --RESTORE_HARDCLIPS false \
	    --RESTORE_ORIGINAL_QUALITIES true \
	    --VALIDATION_STRINGENCY SILENT

	gatk SamToFastq -I "${out_prefix}.mito.ubam" -FASTQ /dev/stdout -INTERLEAVE true \
	    | bwa mem -p -t ${task.cpus} "$ref_fa" /dev/stdin \
	    | samtools sort -@ ${task.cpus} -o "${out_prefix}.mito.aligned.bam" -

	gatk MergeBamAlignment \
	    --UNMAPPED_BAM "${out_prefix}.mito.ubam" \
	    --ALIGNED_BAM "${out_prefix}.mito.aligned.bam" \
	    --REFERENCE_SEQUENCE "$ref_fa" \
	    --OUTPUT "${out_prefix}.mito.merged.bam" \
	    --CREATE_INDEX true \
	    --VALIDATION_STRINGENCY SILENT \
	    --ATTRIBUTES_TO_RETAIN X0

	gatk MarkDuplicates \
	    -I "${out_prefix}.mito.merged.bam" \
	    -O "${out_prefix}.mito.prepped.bam" \
	    -M "${out_prefix}.mito.dup_metrics.txt" \
	    --CREATE_INDEX true \
	    --VALIDATION_STRINGENCY SILENT
	"""
}
