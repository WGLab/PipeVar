
// Phase SNV/SV evidence and aggregate ranked evidence into final prioritized VCF.
process longphase {
	container ='beoungl/docker_test:longphase_0.2.9'

	input:
	tuple path(bam_path), path(index)
	path annovar_vcf
	path annovar_sv_vcf
	path sv_pathogenic
	path snv_rankscore
	path snv_pathogenic
	path snv_rankvar
	path hpo_path
	val out_prefix
	tuple path(ref_fa), path(fa_index)

	output:
	path "${out_prefix}.prio.vcf"

	
	script:


		"""
		/longphase_linux-x64 phase -s $annovar_vcf --sv-file=$annovar_sv_vcf -t 8 -o ${out_prefix}_phased --ont -b $bam_path -r $ref_fa

		/longphase_linux-x64 haplotag -r $ref_fa -s ${out_prefix}_phased.vcf --sv-file ${out_prefix}_phased_SV.vcf -b $bam_path -t 4 -o ${out_prefix}_haplotag

	        bash /clinvar_vcf_and_txt.sh $snv_pathogenic ${out_prefix}_phased.vcf $out_prefix

	        bash /phenosv_vcf_and_tsv.sh $sv_pathogenic ${out_prefix}_phased_SV.vcf $out_prefix

	        bash /rankscore_vcf_and_txt.sh $snv_rankscore ${out_prefix}_phased.vcf $out_prefix

	        bash /rankvar_vcf_and_tsv.sh $snv_rankvar ${out_prefix}_phased.vcf $out_prefix

	python3 /assign_dom_or_rec.py ${out_prefix}.clinvar.vcf ${out_prefix}.phenosv.vcf ${out_prefix}.rankscore.vcf ${out_prefix}.rankvar.vcf $hpo_path OMIM ${out_prefix}.assigned.vcf

	python3 /prio_gene_only.py ${out_prefix}.assigned.vcf ${out_prefix}.prio.vcf




	"""

}
