
// Batch long-read SNP/indel calling with Clair3.
process multi_clair3 {
	container ='hkubal/clair3:v1.2.0'


        input:
	tuple val(out_prefix), path(bam), path(index_file)
	tuple path(ref_fa), path(fa_index)
	val bed_file

	output:
	tuple val(out_prefix), path("${out_prefix}.clair3.vcf.gz")

	script:
	def args   = task.ext.args ?: ''

	"""

	curr_dir=\$PWD

        source /opt/conda/etc/profile.d/conda.sh
	conda activate clair3


	if [ $bed_file != "null" ]; then
		run_clair3.sh --bam_fn=$bam --ref_fn=$ref_fa --threads=4 $args --output=${out_prefix}_clair3 --bed_fn=$bed_file
	else 
		run_clair3.sh --bam_fn=$bam --ref_fn=$ref_fa --threads=4 $args --output=${out_prefix}_clair3
	fi
	
	mv ${out_prefix}_clair3/merge_output.vcf.gz ${out_prefix}_clair3/${out_prefix}.clair3.vcf.gz

	mv ${out_prefix}_clair3/${out_prefix}.clair3.vcf.gz \$curr_dir


	"""


}



