
// Batch phenotype-aware SV scoring with PhenoSV.
process multi_phenosv {
	container ='beoungl/docker_test:phenosv_0.2'


	input:
	tuple val(out_prefix), path(canonical_bed), path(phenosv_bed), path(phenosv_bedpe), path(members_tsv), path(hpo)


	output:
	tuple val(out_prefix), path("${out_prefix}.phenosv.filtered.tsv")

	script:
	def args   = task.ext.args ?: ''
	def min_score = params.phenosv_score ?: '0.50'
	def light_mode = args.contains('PhenoSV-light')

	"""

	source /conda/etc/profile.d/conda.sh
	conda activate phenosv

	HPO_STRING=\$(paste -sd, $hpo)

	curr_dir=\$PWD

	phenosv_dir="\${curr_dir}/${out_prefix}_phenosv"

	mkdir -p \$phenosv_dir

	printf 'SV_ID\tCHROM\tSTART\tEND\tSVTYPE\tPATHOGENICITY\tPHEN2GENE\tPHENOSV_SCORE\tPHENOSV_TYPE\n' > ${out_prefix}.phenosv.simple.tsv
	printf 'SV_ID\tCHROM\tSTART\tEND\tSVTYPE\tPATHOGENICITY\tPHEN2GENE\tPHENOSV_SCORE\tPHENOSV_TYPE\n' > ${out_prefix}.phenosv.bnd.tsv

	simple_rows=\$(awk 'NF {count++} END {print count+0}' $phenosv_bed)
	if [[ \$simple_rows -gt 0 ]]; then
	    mkdir -p \$phenosv_dir/simple
	    python3 /opt/PhenoSV/phenosv/model/phenosv.py --sv_file $phenosv_bed $args --target_folder \$phenosv_dir/simple --target_file_name \$phenosv_dir/simple/phenosv_out --HPO "\$HPO_STRING"
	    python3 /normalize_phenosv_events.py \
	        \$phenosv_dir/simple/phenosv_out.csv \
	        $canonical_bed \
	        ${out_prefix}.phenosv.simple.tsv \
	        --members $members_tsv \
	        --min-score "$min_score"
	fi

	bedpe_rows=\$(awk 'NF {count++} END {print count+0}' $phenosv_bedpe)
	if [[ \$bedpe_rows -gt 0 ]]; then
	    if [[ "${light_mode}" == "true" ]]; then
	        echo "WARNING: PhenoSV-light has reduced accuracy for translocations; scoring BEDPE input as requested." >&2
	    fi
	    mkdir -p \$phenosv_dir/bnd
	    python3 /opt/PhenoSV/phenosv/model/phenosv.py --sv_file $phenosv_bedpe $args --target_folder \$phenosv_dir/bnd --target_file_name \$phenosv_dir/bnd/phenosv_out --HPO "\$HPO_STRING"
	    python3 /normalize_phenosv_events.py \
	        \$phenosv_dir/bnd/phenosv_out.csv \
	        $phenosv_bedpe \
	        ${out_prefix}.phenosv.bnd.tsv \
	        --members $members_tsv \
	        --min-score "$min_score"
	fi

	head -n 1 ${out_prefix}.phenosv.simple.tsv > ${out_prefix}.phenosv.filtered.tsv
	tail -n +2 ${out_prefix}.phenosv.simple.tsv >> ${out_prefix}.phenosv.filtered.tsv
	tail -n +2 ${out_prefix}.phenosv.bnd.tsv >> ${out_prefix}.phenosv.filtered.tsv


	


	"""

}	
