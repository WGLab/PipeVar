// Call mitochondrial SNP/indel variants from long-read alignments using Clair3.
process mito_clair3 {
	container ='hkubal/clair3:v1.2.0'

	input:
	tuple path(bam), path(index)
	val out_prefix
	tuple path(ref_fa), path(fa_index)
	val mito_contig

	output:
	path "${out_prefix}.mito.clair3.raw.vcf.gz"

	script:
	def args = task.ext.args ?: ''

	"""
	source /opt/conda/etc/profile.d/conda.sh
	conda activate clair3

	MITO_CONTIG="${mito_contig}"
	if ! samtools idxstats "$bam" | cut -f1 | grep -qx "\$MITO_CONTIG"; then
	    if samtools idxstats "$bam" | cut -f1 | grep -qx "MT"; then
	        MITO_CONTIG="MT"
	    elif samtools idxstats "$bam" | cut -f1 | grep -qx "chrM"; then
	        MITO_CONTIG="chrM"
	    elif samtools idxstats "$bam" | cut -f1 | grep -qx "chrMT"; then
	        MITO_CONTIG="chrMT"
	    elif samtools idxstats "$bam" | cut -f1 | grep -qx "M"; then
	        MITO_CONTIG="M"
	    else
	        echo "Unable to find mitochondrial contig in $bam" >&2
	        exit 1
	    fi
	fi

	if ! grep -Eq "^>\${MITO_CONTIG}([[:space:]]|\$)" "$ref_fa"; then
	    echo "Mitochondrial contig \${MITO_CONTIG} not found in $ref_fa" >&2
	    exit 1
	fi

	run_clair3.sh \\
	    --bam_fn="$bam" \\
	    --ref_fn="$ref_fa" \\
	    --threads="${task.cpus}" \\
	    $args \\
	    --output="${out_prefix}_mito_clair3" \\
	    --ctg_name="\$MITO_CONTIG" \\
	    --haploid_sensitive \\
	    --snp_min_af=0.01 \\
	    --indel_min_af=0.05 \\
	    --var_pct_full=0.1 \\
	    --enable_variant_calling_at_sequence_head_and_tail

	mv "${out_prefix}_mito_clair3/merge_output.vcf.gz" "${out_prefix}.mito.clair3.raw.vcf.gz"
	"""
}
