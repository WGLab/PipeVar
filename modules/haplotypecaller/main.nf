
process haplotypecaller {
	container='broadinstitute/gatk:4.5.0.0'

        input:
        tuple path(bam), path(index_file)
        val out_prefix
	tuple path(ref_fa), path(fa_index), path(dict_index)
	val bed_file

	output:
	path "${out_prefix}.recal.vcf.gz"

	script:

	"""


	if [ $bed_file != "null" ]; then
		gatk HaplotypeCaller --input $bam --output ${out_prefix}.vcf.gz --reference $ref_fa -L $bed_file
	else
		gatk HaplotypeCaller --input $bam --output ${out_prefix}.vcf.gz --reference $ref_fa
	fi

	#Download all necessary file here, and delete them later on
	
	mkdir -p ./gatk_files	

	wget -O "./gatk_files/1000G_omni2.5.hg38.vcf.gz" https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/1000G_omni2.5.hg38.vcf.gz
	wget -O "./gatk_files/1000G_omni2.5.hg38.vcf.gz.tbi" https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/1000G_omni2.5.hg38.vcf.gz.tbi
	wget -O "./gatk_files/1000G_phase1.snps.high_confidence.hg38.vcf.gz" https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz
	wget -O "./gatk_files/1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi" https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi
	wget -O "./gatk_files/Homo_sapiens_assembly38.dbsnp138.vcf" https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf
	bgzip -f ./gatk_files/Homo_sapiens_assembly38.dbsnp138.vcf
	tabix -f ./gatk_files/Homo_sapiens_assembly38.dbsnp138.vcf.gz
	wget -O "./gatk_files/hapmap_3.3.hg38.vcf.gz" https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/hapmap_3.3.hg38.vcf.gz
	wget -O "./gatk_files/hapmap_3.3.hg38.vcf.gz.tbi" https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/hapmap_3.3.hg38.vcf.gz.tbi
	

	gatk VariantRecalibrator -R $ref_fa -V ${out_prefix}.vcf.gz --resource:hapmap,known=false,training=true,truth=true,prior=15.0 ./gatk_files/hapmap_3.3.hg38.vcf.gz --resource:omni,known=false,training=true,truth=false,prior=12.0 ./gatk_files/1000G_omni2.5.hg38.vcf.gz --resource:1000G,known=false,training=true,truth=false,prior=10.0 ./gatk_files/1000G_phase1.snps.high_confidence.hg38.vcf.gz --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ./gatk_files/Homo_sapiens_assembly38.dbsnp138.vcf.gz -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -mode SNP -O ${out_prefix}.recal --tranches-file ${out_prefix}.tranches

	 gatk ApplyVQSR -R $ref_fa -V ${out_prefix}.vcf.gz -O ${out_prefix}.recal.vcf.gz --ts-filter-level 99.0 --tranches-file ${out_prefix}.tranches --recal-file ${out_prefix}.recal -mode SNP
 
 
	rm -r ./gatk_files
	

	"""


}



