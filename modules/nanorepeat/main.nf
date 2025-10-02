
process NanoRepeat {
        container ='beoungl/docker_test:nanorepeat_beta'


        input:
        path bam
	val input_directory
        val out_prefix
	path ref_fa
	val ref_fa_directory
	path output_directory
	val output_directory_full




	script:
	def args  = task.ext.args ?: ''	

	"""

	python3 /usr/local/lib/python3.10/dist-packages/NanoRepeat/nanoRepeat.py -i $input_directory/$bam -t bam -d $args -r $ref_fa_directory/$ref_fa -b /Nanorepeat_bed/nanorepeat.input.bed -o $output_directory_full/${out_prefix}_nanoRepeat_output

	sh /Nanorepeat_bed/compare_nanorepeat.sh $output_directory_full/${out_prefix}_nanoRepeat_output.NanoRepeat_output.tsv > $output_directory_full/${out_prefix}_nanorepeat_result.tsv





	"""


}



