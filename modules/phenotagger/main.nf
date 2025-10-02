
process phenotagger {
	container ='beoungl/docker_test:phenotagger'


        input:
        path phenotagger_input
	val phenotagger_directory
	path output_directory
	val output_directory_full
	val output_prefix


	output:
	val "$output_directory_full/${output_prefix}_phenotagger_patient_hpo.txt"

	"""

	python3 /PhenoTagger_py/generate_input.py -i $phenotagger_directory/$phenotagger_input -t $output_directory_full/${output_prefix}_phenotagger/phenotagger_temp -o $output_directory_full/${output_prefix}_phenotagger

 	cd /opt/PhenoTagger/src

	python3 PhenoTagger_tagging.py -i $output_directory_full/${output_prefix}_phenotagger/phenotagger_temp/phenotagger_input/input0/ -o $output_directory_full/${output_prefix}_phenotagger/phenotagger_output/output0/

	awk -F'\t' '\$6 ~ /HP:/ {print \$6}' $output_directory_full/${output_prefix}_phenotagger/phenotagger_output/output0/patient_note.PubTator > $output_directory_full/${output_prefix}_phenotagger_patient_hpo.txt



	"""


}



