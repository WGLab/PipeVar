
process PhenoSV {
        container ='beoungl/docker_test:phenosv'


	input:
	path bed
	val out_prefix
	path hpo
	val hpo_directory
	path input_directory_path
	val input_directory


	output:
	val "${out_prefix}.phenosv.filtered.tsv"

	script:
	def args   = task.ext.args ?: ''

	"""

	source /conda/etc/profile.d/conda.sh
	conda activate phenosv

	HPO_STRING=\$(paste -sd, $hpo_directory/$hpo)

	mkdir -p $input_directory/${out_prefix}_phenosv


	python3 /opt/PhenoSV/phenosv/model/phenosv.py --sv_file $bed $args --target_folder $input_directory/${out_prefix}_phenosv --target_file_name  $input_directory/${out_prefix}_phenosv/phenosv_out --HPO "\$HPO_STRING"

	awk -F',' '\$6 > 0.5' $input_directory/${out_prefix}_phenosv/phenosv_out.csv | awk -F',' '\$2 == "SV" ' | awk -F',' '{print \$7"\t"\$0}' | sort -k1,1 > $input_directory/phenosv_top.join.tsv

	sort -k4,4 $input_directory/${out_prefix}.bed > $input_directory/${out_prefix}.sorted.bed


	join -t\$'\t' -1 1 -2 4 $input_directory/phenosv_top.join.tsv $input_directory/${out_prefix}.sorted.bed | awk -F'\t' 'BEGIN{OFS="\t"} { id=\$1; \$1=""; sub(/^\t/, ""); print \$0, id}' > $input_directory/${out_prefix}.phenosv.filtered.tsv

	rm $input_directory/phenosv_top.join.tsv
	rm $input_directory/${out_prefix}.sorted.bed


	


	"""

}	
