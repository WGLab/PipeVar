
process multi_haplotypecaller {
	container='broadinstitute/gatk:4.5.0.0'

        input:
	tuple val(out_prefix), path(bam), path(index_file), path(omni_vcf), path(omni_vcf_tabix), path(phase_snp_vcf), path(phase_snp_vcf_tbi), path(dbsnp_vcf), path(dbsnp_vcf_tbi), path(hapmap_vcf), path(hapmap_vcf_tbi)
	tuple path(ref_fa), path(fa_index), path(dict_index)
	val bed_file

	output:
	tuple val(out_prefix), path("${out_prefix}.recal.vcf.gz")

	script:

	"""

	if [ $bed_file != "null" ]; then
		gatk HaplotypeCaller --input $bam --output ${out_prefix}.vcf.gz --reference $ref_fa -L $bed_file
	else
		gatk HaplotypeCaller --input $bam --output ${out_prefix}.vcf.gz --reference $ref_fa
	fi

	
	gatk VariantRecalibrator -R $ref_fa -V ${out_prefix}.vcf.gz --resource:hapmap,known=false,training=true,truth=true,prior=15.0 $hapmap_vcf --resource:omni,known=false,training=true,truth=false,prior=12.0 $omni_vcf --resource:1000G,known=false,training=true,truth=false,prior=10.0 $phase_snp_vcf --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $hapmap_vcf -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -mode SNP -O ${out_prefix}.recal --tranches-file ${out_prefix}.tranches

	 gatk ApplyVQSR -R $ref_fa -V ${out_prefix}.vcf.gz -O ${out_prefix}.recal.vcf.gz --ts-filter-level 99.0 --tranches-file ${out_prefix}.tranches --recal-file ${out_prefix}.recal -mode SNP
 
 
	

	"""


}



