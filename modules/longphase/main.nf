
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


	/longphase_linux-x64 phase -s $output_directory_full/$snv_vcf_path --sv-file=$output_directory_full/$sv_vcf_path -t 8 -o $output_directory_full/$out_prefix --ont -b $bam_directory/$bam_path -r $ref_fa_full/$ref_fa

	#Process long phase vcf results here, given the phasing.
	#Based on ACMG guideline.


	#Collect rankscore vcf + rankvar analysis here and match with phasing information.

	#Input would be phased snv file, rankvar result with certain filter.

	#sh snv_phase_rankvar.sh $snv_rankvar ${out_prefix}.vcf

	#Input would be phased snv file, rankscore analysis with certain filter.

	#sh snv_phase_rankscore.sh $snv_rankscore ${out_prefix}.vcf

	#Input would be phased snv file with Clinvar.

	#sh snv_phase_pathogenic.sh $snv_pathogenic ${out_prefix}.vcf

	#Input would be phased sv file with PhenoSV filtered result.

	#sh sv_phase_pathogenic.sh $sv_pathogenic ${out_prefix}_SV.vcf

	#All the output should be just filtered sv/snv based on the result.

	#Combined SNV/SV script that combines all the results here.

	#First, add the ClinVar result into the top.

	#See if you can match with the toerh.
	


	#Goes in the order of
	#Pathogenic variants, known in ClinVar, AD.
	#Pathgoenic variants, known in ClinVar, AR.
	#X-linked inheritance somewhere in between.
	#Pathogenic variants + VUS, AR.
	#VUS, AD.
	#VUS + VUS, AR.
	#VUS, AR.


	"""

}

