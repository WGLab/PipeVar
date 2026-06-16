
// Batch longphase phasing and evidence aggregation into prioritized VCFs.
process multi_longphase {
	container ='beoungl/docker_test:longphase_0.2.29'

	input:
	tuple val(out_prefix), path(snv_rankvar), path(snv_rankscore), path(snv_pathogenic), path(sv_pathogenic), path(sv_vcf_path), path(snv_vcf_path), path(bam_path), path(bam_index), path(hpo_path), val(age_of_onset), val(sex)
	tuple path(ref_fa), path(fa_index)

	val(inheritance_mode)
	val(include_clinvar_report)
	val(allow_unphased_comphet)

	output:
	tuple path("${out_prefix}.prio.vcf"), path("${out_prefix}.prio_gene.vcf"), path("${out_prefix}_haplotag.bam")
	
	script:
	def platform_args = task.ext.args ?: '--ont'
	def min_score = params.phenosv_score ?: '0.50'
	def gene_filter = params.gene ?: ''
	def sv_only = params.prioritize_sv_only ?: 'no'


	"""
	#Final process to organize the results and show them in ACMG guideline format; might take long to develop.

	#Join the bam path, vcfs and output files


		/longphase_linux-x64 phase -s $snv_vcf_path --sv-file=$sv_vcf_path -t ${task.cpus} -o ${out_prefix}_phased $platform_args -b $bam_path -r $ref_fa

		/longphase_linux-x64 haplotag -r $ref_fa -s ${out_prefix}_phased.vcf --sv-file ${out_prefix}_phased_SV.vcf -b $bam_path -t ${task.cpus} -o ${out_prefix}_haplotag

	        bash /phenosv_vcf_and_tsv.sh $sv_pathogenic ${out_prefix}_phased_SV.vcf $out_prefix
		mv ${out_prefix}.phenosv.vcf ${out_prefix}.phenosv.unfiltered.vcf
		python3 /filter_phenosv_vcf.py ${out_prefix}.phenosv.unfiltered.vcf ${out_prefix}.phenosv.vcf --min-score $min_score --genes "$gene_filter"

	if [[ "$sv_only" == "yes" ]]; then
	    python3 /assign_dom_or_rec_sv_only.py ${out_prefix}.phenosv.vcf $hpo_path $age_of_onset $inheritance_mode ${out_prefix}.assigned.vcf --phenosv-score $min_score --genes "$gene_filter"
	else
	    bash /clinvar_vcf_and_txt.sh $snv_pathogenic ${out_prefix}_phased.vcf $out_prefix
	    bash /rankscore_vcf_and_txt.sh $snv_rankscore ${out_prefix}_phased.vcf $out_prefix
	    bash /rankvar_vcf_and_tsv.sh $snv_rankvar ${out_prefix}_phased.vcf $out_prefix
	    python3 /assign_dom_or_rec.py ${out_prefix}.clinvar.vcf ${out_prefix}.phenosv.vcf ${out_prefix}.rankscore.vcf ${out_prefix}.rankvar.vcf $hpo_path $age_of_onset $inheritance_mode ${out_prefix}.assigned.vcf --phenosv-score $min_score --genes "$gene_filter"
	fi

	python3 /prio_gene_only.py ${out_prefix}.assigned.vcf ${out_prefix}.prio_gene.vcf gene --include-clinvar $include_clinvar_report --allow-unphased-comphet $allow_unphased_comphet --genes "$gene_filter" --sex $sex
	python3 /prio_gene_only.py ${out_prefix}.assigned.vcf ${out_prefix}.prio.vcf variant --include-clinvar $include_clinvar_report --allow-unphased-comphet $allow_unphased_comphet --genes "$gene_filter" --sex $sex


	"""

}
