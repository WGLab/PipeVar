
// Batch long-read SNP/indel calling with Clair3.
process multi_clair3 {
	container ='hkubal/clair3:v1.2.0'


        input:
	tuple val(out_prefix), path(bam), path(index_file), path(bed_file)
	tuple path(ref_fa), path(fa_index)

	output:
	tuple val(out_prefix), path("${out_prefix}.clair3.vcf.gz")

	script:
	def args   = task.ext.args ?: ''
	def regions = bed_file ? "--bed_fn=${bed_file}" : ""

	"""

	curr_dir=\$PWD

        source /opt/conda/etc/profile.d/conda.sh
	conda activate clair3


	run_clair3.sh --bam_fn=$bam --ref_fn=$ref_fa --threads=${task.cpus} $args --output=${out_prefix}_clair3 $regions
	
	mv ${out_prefix}_clair3/merge_output.vcf.gz ${out_prefix}_clair3/${out_prefix}.clair3.vcf.gz

	mv ${out_prefix}_clair3/${out_prefix}.clair3.vcf.gz \$curr_dir


	"""


}

