
// Score and prioritize structural variants with phenotype-aware PhenoSV model.
process PhenoSV {
        container ='beoungl/docker_test:phenosv'


	input:
	path bed
	val out_prefix
	path hpo


	output:
	path "${out_prefix}.phenosv.filtered.tsv"

	script:
	def args   = task.ext.args ?: ''

	"""

	source /conda/etc/profile.d/conda.sh
	conda activate phenosv
	
	curr_dir=\$PWD

	phenosv_dir="\${curr_dir}/${out_prefix}_phenosv"

	HPO_STRING=\$(paste -sd, $hpo)

	mkdir -p \$phenosv_dir

	# If upstream SV filtering produced no records, emit an empty output and continue.
	# Header/comment-only BED should be treated as empty input.
	data_rows=\$(awk 'NF && \$0 !~ /^#/' $bed | wc -l)
	if [[ ! -s $bed ]] || [[ \$data_rows -eq 0 ]]; then
	    : > ${out_prefix}.phenosv.filtered.tsv
	    exit 0
	fi


	python3 /opt/PhenoSV/phenosv/model/phenosv.py --sv_file $bed $args --target_folder ${out_prefix}_phenosv --target_file_name  \$phenosv_dir/phenosv_out --HPO "\$HPO_STRING"

	awk -F',' '\$6 > 0.5 && \$7 != ""' \$phenosv_dir/phenosv_out.csv | awk -F',' '{print \$7"\t"\$0}' | sort -k1,1 > ${out_prefix}_phenosv_top.join.tsv || true

	sort -k4,4 $bed > ${out_prefix}.sorted.bed

	if [[ ! -s ${out_prefix}_phenosv_top.join.tsv ]]; then
	    : > ${out_prefix}.phenosv.filtered.tsv
	else

	join -t\$'\t' -1 1 -2 4 ${out_prefix}_phenosv_top.join.tsv ${out_prefix}.sorted.bed | awk -F'\t' 'BEGIN{OFS="\t"} { id=\$1; \$1=""; sub(/^\t/, ""); print \$0, id}' | sed 's/,/\t/g' | awk -F'\t' 'BEGIN{OFS="\t"} {print \$8,\$9,\$10,\$12,\$3,\$4,\$5,\$6}' > ${out_prefix}.phenosv.filtered.tsv || true

	# `join` returns non-zero when there are no overlaps; keep pipeline alive with empty output.
	if [[ ! -s ${out_prefix}.phenosv.filtered.tsv ]]; then
	    : > ${out_prefix}.phenosv.filtered.tsv
	fi
	fi

	rm -f ${out_prefix}_phenosv_top.join.tsv
	rm -f ${out_prefix}.sorted.bed


	


	"""

}	
