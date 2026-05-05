
// Call SNP/indel variants from long-read alignments using NanoCaller.
process nanocaller {
        container ='beoungl/docker_test:nanocaller'


        input:
	tuple path(bam), path(index)
        val out_prefix
	tuple path(ref_fa), path(fa_index)
	val bed_file

	output:
	path "${out_prefix}.nanocaller.vcf.gz"
	
	script:
	def args   = task.ext.args ?: ''


	"""
	source /conda/etc/profile.d/conda.sh
        conda activate nanocaller

	if [ $bed_file != "null" ]; then
		NanoCaller --bam $bam --ref $ref_fa --cpu ${task.cpus} --output $out_prefix $args --bed $bed_file
	else
		NanoCaller --bam $bam --ref $ref_fa --cpu ${task.cpus} --output $out_prefix $args
	fi

	mv $out_prefix/variant_calls.vcf.gz $out_prefix/${out_prefix}.nanocaller.vcf.gz

	curr_dir=\$PWD

	mv $out_prefix/${out_prefix}.nanocaller.vcf.gz \$curr_dir


	"""


}


