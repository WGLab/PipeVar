
// Batch extraction of HPO terms from clinical notes using PhenoTagger.
process multi_phenotagger {
	container ='beoungl/docker_test:phenotagger'


        input:
	tuple val(out_prefix),  path(phenotagger_input)


	output:
	tuple val(out_prefix), path("${out_prefix}_phenotagger_patient_hpo.txt")

	"""
	export HOME=\$PWD


	curr_dir=\$PWD

	temp_dir="\${curr_dir}/${out_prefix}_phenotagger"

	mv $phenotagger_input ${out_prefix}_phenotagger_input.txt	

	python3 /PhenoTagger_py/generate_input.py \
        -i ${out_prefix}_phenotagger_input.txt \
        -t \${temp_dir}/phenotagger_temp \
        -o \${temp_dir}




 	cd /opt/PhenoTagger/src

	
	python3 PhenoTagger_tagging.py \
        -i \${temp_dir}/phenotagger_temp/phenotagger_input/input0/ \
        -o \${temp_dir}/phenotagger_output/output0/
	
	awk -F'\t' '\$6 ~ /HP:/ {print \$6}' \
        \${temp_dir}/phenotagger_output/output0/${out_prefix}_phenotagger_input.PubTator > \${curr_dir}/${out_prefix}_phenotagger_patient_hpo.txt


	"""


}



