
// Batch repeat expansion calling from long-read BAM inputs with NanoRepeat.
process multi_nanorepeat {
        container ='beoungl/docker_test:nanorepeat_0.1'


        input:
	tuple val(out_prefix), path(bam), path(index_file)
	tuple path(ref_fa), path(fa_index)

	output:
	path "${out_prefix}_nanorepeat_result.tsv"
	


	script:
	def args  = task.ext.args ?: ''	

	"""

	nanoRepeat.py -i $bam -t bam -d $args -r $ref_fa -b /Nanorepeat_bed/nanorepeat.input.bed -o ${out_prefix}_nanoRepeat_output

	sh /Nanorepeat_bed/compare_nanorepeat.sh ${out_prefix}_nanoRepeat_output.NanoRepeat_output.tsv > ${out_prefix}_nanorepeat_result.tsv





	"""


}



