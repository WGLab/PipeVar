
process longphase {
	container='beoungl/docker_test:longphase'

	input:
	path bam_path
	val bam_directory
	path snv_vcf_path
	path sv_vcf_path
	path output_directory
	val output_directory_full
	val sv_pathogenic
	val snv_rankscore
	val snv_pathogenic
	val snv_rankvar
	val out_prefix
	path ref_fa
	val ref_fa_full



	
	script:


	"""
	#Final process to organize the results and show them in ACMG guideline format; might take long to develop.


	/longphase_linux-x64 phase -s $output_directory_full/$snv_vcf_path --sv-file=$output_directory_full/${sv_vcf_path}.gz -t 8 -o $output_directory_full/$out_prefix --ont -b $bam_directory/$bam_path -r $ref_fa_full/$ref_fa

	/longphase_linux-x64 haplotag -r $ref_fa_full/$ref_fa -s $output_directory_full/${out_prefix}.vcf --sv-file $output_directory_full/${out_prefix}_SV.vcf -b $bam_directory/$bam_path -t 4 -o $output_directory_full/${out_prefix}_haplotag



	"""

}

