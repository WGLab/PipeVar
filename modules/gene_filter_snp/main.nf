
process gene_filter_snp {
	container ='/mnt/isilon/wang_lab/beounglee/projects/pipeline/config_testing/modules_images/clair3.sif'


        input:
	path vcf
	val input_directory
	val gene_infomration
        val out_prefix
	path output_directory
	val output_directory_full

	output:


	"""

	grep '$gene_information' $vcf > 	


	"""


}



