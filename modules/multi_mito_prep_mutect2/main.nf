process multi_mito_prep_mutect2 {
	container = 'beoungl/docker_test:mito_mutect2_0.1'

	input:
	tuple val(out_prefix), path(bam), path(index_file)
	tuple path(ref_fa), path(fa_index), path(dict_index), path(bwa_amb), path(bwa_ann), path(bwa_bwt), path(bwa_pac), path(bwa_sa)
	val mito_contig

	output:
	tuple val(out_prefix), path("${out_prefix}.mito.prepped.bam"), path("${out_prefix}.mito.prepped.bam.bai"), path("${out_prefix}.mito.dup_metrics.txt")

	script:
	"""
	MITO_CONTIG="${mito_contig}"
	CONTIGS="\$(samtools view -H "$bam" | awk 'BEGIN {FS="\\t"} /^@SQ/ {for (i=1; i<=NF; i++) if (\$i ~ /^SN:/) {sub(/^SN:/, "", \$i); print \$i}}')"

	if ! printf '%s\\n' "\$CONTIGS" | grep -qx "\$MITO_CONTIG"; then
	    if printf '%s\\n' "\$CONTIGS" | grep -qx "MT"; then
	        MITO_CONTIG="MT"
	    elif printf '%s\\n' "\$CONTIGS" | grep -qx "chrM"; then
	        MITO_CONTIG="chrM"
	    elif printf '%s\\n' "\$CONTIGS" | grep -qx "M"; then
	        MITO_CONTIG="M"
	    elif printf '%s\\n' "\$CONTIGS" | grep -qx "chrMT"; then
	        MITO_CONTIG="chrMT"
	    else
	        echo "Unable to find mitochondrial contig in $bam" >&2
	        exit 1
	    fi
	fi

	samtools view -@ ${task.cpus} -b -T "$ref_fa" -o "${out_prefix}.mito.subset.unsorted.bam" "$bam" "\$MITO_CONTIG"
	samtools sort -@ ${task.cpus} -o "${out_prefix}.mito.subset.bam" "${out_prefix}.mito.subset.unsorted.bam"
	samtools index -@ ${task.cpus} "${out_prefix}.mito.subset.bam"
	INPUT_BAM="${out_prefix}.mito.subset.bam"

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
	fi

	gatk RevertSam \
	    -I "\$INPUT_BAM" \
	    -O "${out_prefix}.mito.ubam" \
	    --REMOVE_ALIGNMENT_INFORMATION true \
	    --RESTORE_HARDCLIPS false \
	    --RESTORE_ORIGINAL_QUALITIES true \
	    --SANITIZE true \
	    --MAX_DISCARD_FRACTION 0.05 \
	    --SORT_ORDER queryname \
	    --VALIDATION_STRINGENCY LENIENT

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

	mv "${out_prefix}.mito.prepped.bai" "${out_prefix}.mito.prepped.bam.bai"
	rm -f \
	    "${out_prefix}.mito.subset.unsorted.bam" \
	    "${out_prefix}.mito.subset.bam" \
	    "${out_prefix}.mito.subset.bam.bai" \
	    "${out_prefix}.mito.rg.bam" \
	    "${out_prefix}.mito.rg.bai" \
	    "${out_prefix}.mito.rg.bam.bai" \
	    "${out_prefix}.mito.ubam" \
	    "${out_prefix}.mito.aligned.bam" \
	    "${out_prefix}.mito.merged.bam" \
	    "${out_prefix}.mito.merged.bai" \
	    "${out_prefix}.mito.merged.bam.bai"
	"""
}
