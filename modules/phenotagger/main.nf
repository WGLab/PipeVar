
// Extract HPO terms from unstructured clinical notes using PhenoTagger.
process phenotagger {
	
	container ='beoungl/docker_test:phenotagger'


        input:
        path phenotagger_input
	val output_prefix


	output:
	path "${output_prefix}_phenotagger_patient_hpo.txt"

	"""
	#Need to do this because it seems like tensorflow thinks we need to be at home directory for some reason
	export HOME=\$PWD

	curr_dir=\$PWD

	temp_dir="\${curr_dir}/${output_prefix}_phenotagger"

	normalized_input=${output_prefix}_phenotagger_input.normalized.tmp
	awk '{
		gsub(/\r/, " ", \$0)
		for (i = 1; i <= NF; i++) {
			printf "%s%s", separator, \$i
			separator = " "
		}
	}
	END {
		if (separator != "") {
			print ""
		}
	}' $phenotagger_input > \${normalized_input}
	mv \${normalized_input} ${output_prefix}_phenotagger_input.txt


	python3 /PhenoTagger_py/generate_input.py \
   	-i ${output_prefix}_phenotagger_input.txt \
        -t \${temp_dir}/phenotagger_temp \
        -o \${temp_dir}


	cd /opt/PhenoTagger/src


	python3 PhenoTagger_tagging.py \
        -i \${temp_dir}/phenotagger_temp/phenotagger_input/input0/ \
        -o \${temp_dir}/phenotagger_output/output0/

	awk -F'\t' '\$6 ~ /HP:/ {print \$6}' \
        \${temp_dir}/phenotagger_output/output0/${output_prefix}_phenotagger_input.PubTator > \${curr_dir}/${output_prefix}_phenotagger_patient_hpo.txt

	"""


}


