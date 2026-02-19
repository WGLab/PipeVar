
process multi_longphase {
	container='beoungl/docker_test:longphase_0.2.8'

	input:
	tuple val(out_prefix), path(snv_rankvar), path(snv_rankscore), path(snv_pathogenic), path(sv_pathogenic), path(sv_vcf_path), path(snv_vcf_path), path(bam_path), path(bam_index)
	tuple path(ref_fa), path(fa_index)


	output:
	tuple path("${out_prefix}.prio.vcf"), path("${out_prefix}_haplotag.bam")
	
	script:


	"""
	#Final process to organize the results and show them in ACMG guideline format; might take long to develop.

	#Join the bam path, vcfs and output files


	/longphase_linux-x64 phase -s $snv_vcf_path --sv-file=$sv_vcf_path -t 8 -o ${out_prefix}_phased --ont -b $bam_path -r $ref_fa

	/longphase_linux-x64 haplotag -r $ref_fa -s ${out_prefix}_phased.vcf --sv-file ${out_prefix}_phased_SV.vcf -b $bam_path -t 4 -o ${out_prefix}_haplotag


        bash /clinvar_vcf_and_txt.sh $snv_pathogenic ${out_prefix}_phased.vcf $out_prefix

        bash /phenosv_vcf_and_tsv.sh $sv_pathogenic ${out_prefix}_phased_SV.vcf $out_prefix

        bash /rankscore_vcf_and_txt.sh $snv_rankscore ${out_prefix}_phased.vcf $out_prefix

        bash /rankvar_vcf_and_tsv.sh $snv_rankvar ${out_prefix}_phased.vcf $out_prefix

	python3 /assign_dom_or_rec.py ${out_prefix}.clinvar.vcf ${out_prefix}.phenosv.vcf ${out_prefix}.rankscore.vcf ${out_prefix}.rankvar.vcf OMIM ${out_prefix}.assigned.vcf

	python3 /prio_gene_only.py ${out_prefix}.assigned.vcf ${out_prefix}.prio.vcf


	"""

}

