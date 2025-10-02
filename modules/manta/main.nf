
process Manta {
	container ='beoungl/docker_test:manta'


        input:
        path bam
        val input_directory
	val out_prefix
	path ref_fa
	val ref_fa_directory
	path output_directory
	val output_directory_full


	output:
	val "$output_directory_full/${out_prefix}_manta.vcf"

	script:

	"""

	/manta/bin/configManta.py --bam=$input_directory/$bam --referenceFasta=$ref_fa_directory/$ref_fa --runDir $output_directory_full/${out_prefix}_manta
	
	$output_directory_full/${out_prefix}_manta/runWorkflow.py -j 4 -g 64 -m local
	
	gunzip $output_directory_full/${out_prefix}_manta/results/variants/diploidSV.vcf.gz

	mv $output_directory_full/${out_prefix}_manta/results/variants/diploidSV.vcf $output_directory_full/${out_prefix}_manta.vcf

	"""


}



