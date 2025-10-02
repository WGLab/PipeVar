
process Phen2gene {
	container ='beoungl/docker_test:rankvar'


        input:
        path hpo
	val hpo_directory
        val out_prefix
	path input_directory_path
	val input_directory

	output:
	val "$input_directory/${out_prefix}_phen2gene/output_file.associated_gene_list"

	script:

	"""
	source /conda/etc/profile.d/conda.sh
        conda activate phen2gene

	mkdir -p $input_directory/${out_prefix}_phen2gene

	python3 /opt/Phen2Gene/phen2gene.py -f $hpo_directory/$hpo -out $input_directory/${out_prefix}_phen2gene

	"""


}



