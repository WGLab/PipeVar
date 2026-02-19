
// Integrate NGS SNV/SV evidence and produce a prioritized VCF.
process ngs_prio {
	container ='beoungl/docker_test:longphase_0.2.8'

	input:
	val(out_prefix)
	path(snv_rankvar)
	path(snv_rankscore)
	path(snv_pathogenic)
	path(sv_pathogenic)
	path(sv_vcf_path)
	path(snv_vcf_path)

	output:
	path "${out_prefix}.prio.vcf"

	
	script:


	"""

	#No bam file processing, so only sv + snp processing.

        bash /clinvar_vcf_and_txt.sh $snv_pathogenic $snv_vcf_path $out_prefix

        bash /phenosv_vcf_and_tsv.sh $sv_pathogenic $sv_vcf_path $out_prefix

        bash /rankscore_vcf_and_txt.sh $snv_rankscore $snv_vcf_path $out_prefix

        bash /rankvar_vcf_and_tsv.sh $snv_rankvar $snv_vcf_path $out_prefix

	python3 /assign_dom_or_rec.py ${out_prefix}.clinvar.vcf ${out_prefix}.phenosv.vcf ${out_prefix}.rankscore.vcf ${out_prefix}.rankvar.vcf OMIM ${out_prefix}.assigned.vcf

	python3 /prio_gene_only.py ${out_prefix}.assigned.vcf ${out_prefix}.prio.vcf

	"""

}

