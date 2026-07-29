
// Batch long-read SNP/indel calling with NanoCaller.
process multi_nanocaller {
        container ='beoungl/docker_test:nanocaller'


        input:
	tuple val(out_prefix), path(bam), path(index_file), path(bed_file)
	tuple path(ref_fa), path(fa_index)

	output:
	tuple val(out_prefix), path("${out_prefix}.nanocaller.vcf.gz")
	
	script:
	def args   = task.ext.args ?: ''
	def regions = bed_file ? "--bed ${bed_file}" : ""


	"""
	source /conda/etc/profile.d/conda.sh
        conda activate nanocaller

	NanoCaller --bam $bam --ref $ref_fa --cpu ${task.cpus} --output $out_prefix $args $regions


	mv $out_prefix/variant_calls.vcf.gz $out_prefix/${out_prefix}.nanocaller.vcf.gz

	curr_dir=\$PWD

	mv $out_prefix/${out_prefix}.nanocaller.vcf.gz \$curr_dir


	"""


}

