
process multi_phen2gene {
	container ='beoungl/docker_test:rankvar'


        input:
        tuple val(out_prefix), path(note)

	output:
	tuple val(out_prefix), path ("${out_prefix}.phen2gene.txt")

	script:

	"""
	source /conda/etc/profile.d/conda.sh
        conda activate phen2gene

	mkdir -p ${out_prefix}_phen2gene

	python3 /opt/Phen2Gene/phen2gene.py -f $note -out ${out_prefix}_phen2gene

	mv ${out_prefix}_phen2gene/output_file.associated_gene_list ${out_prefix}.phen2gene.txt

	"""


}



