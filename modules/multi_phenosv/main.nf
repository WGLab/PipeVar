
// Batch phenotype-aware SV scoring with PhenoSV.
process multi_phenosv {
        container ='beoungl/docker_test:phenosv'


	input:
	tuple val(out_prefix), path(bed), path(hpo)


	output:
	tuple val(out_prefix), path("${out_prefix}.phenosv.filtered.tsv")

	script:
	def args   = task.ext.args ?: ''

	"""

	source /conda/etc/profile.d/conda.sh
	conda activate phenosv

	HPO_STRING=\$(paste -sd, $hpo)

	curr_dir=\$PWD

	phenosv_dir="\${curr_dir}/${out_prefix}_phenosv"

	mkdir -p \$phenosv_dir


	python3 /opt/PhenoSV/phenosv/model/phenosv.py --sv_file $bed $args --target_folder ${out_prefix}_phenosv --target_file_name  \$phenosv_dir/phenosv_out --HPO "\$HPO_STRING"

	awk -F',' '\$6 > 0.5' \$phenosv_dir/${out_prefix}_phenosv/phenosv_out.csv | awk -F',' '\$2 == "SV" ' | awk -F',' '{print \$7"\t"\$0}' | sort -k1,1 > ${out_prefix}_phenosv_top.join.tsv

	sort -k4,4 ${out_prefix}.bed > ${out_prefix}.sorted.bed


	join -t\$'\t' -1 1 -2 4 ${out_prefix}_phenosv_top.join.tsv ${out_prefix}.sorted.bed | awk -F'\t' 'BEGIN{OFS="\t"} { id=\$1; \$1=""; sub(/^\t/, ""); print \$0, id}' > ${out_prefix}.phenosv.filtered.tsv

	rm ${out_prefix}_phenosv_top.join.tsv
	rm ${out_prefix}.sorted.bed


	


	"""

}	
