
// Batch phenotype-aware SV scoring with PhenoSV.
process multi_phenosv {
        container ='beoungl/docker_test:phenosv_0.1'


	input:
	tuple val(out_prefix), path(bed), path(hpo)


	output:
	tuple val(out_prefix), path("${out_prefix}.phenosv.filtered.tsv")

	script:
	def args   = task.ext.args ?: ''
	def min_score = params.phenosv_score ?: '0.50'

	"""

	source /conda/etc/profile.d/conda.sh
	conda activate phenosv

	HPO_STRING=\$(paste -sd, $hpo)

	curr_dir=\$PWD

	phenosv_dir="\${curr_dir}/${out_prefix}_phenosv"

	mkdir -p \$phenosv_dir

	# Always emit the canonical event TSV header, including for empty SV input.
	data_rows=\$(awk 'NF && \$0 !~ /^#/' $bed | wc -l)
	if [[ ! -s $bed ]] || [[ \$data_rows -eq 0 ]]; then
	    printf 'SV_ID\tCHROM\tSTART\tEND\tSVTYPE\tPATHOGENICITY\tPHEN2GENE\tPHENOSV_SCORE\tPHENOSV_TYPE\n' > ${out_prefix}.phenosv.filtered.tsv
	    exit 0
	fi


	python3 /opt/PhenoSV/phenosv/model/phenosv.py --sv_file $bed $args --target_folder ${out_prefix}_phenosv --target_file_name  \$phenosv_dir/phenosv_out --HPO "\$HPO_STRING"

	python3 /normalize_phenosv_events.py \
	    \$phenosv_dir/phenosv_out.csv \
	    $bed \
	    ${out_prefix}.phenosv.filtered.tsv \
	    --min-score "$min_score"


	


	"""

}	
