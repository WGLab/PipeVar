
process haplotypecaller {
	container='broadinstitute/gatk:4.5.0.0'

        input:
        path bam
	val input_directory
        val out_prefix
	path ref_fa
	val ref_fa_directory
	path output_directory
	val output_directory_full

	output:
	path "$output_directory/${out_prefix}.recal.vcf.gz"

	script:

	"""

	gatk HaplotypeCaller --input $input_directory/$bam --output $output_directory_full/${out_prefix}.vcf.gz --reference $ref_fa_directory/$ref_fa

	#Download all necessary file here, and delete them later on
	
	mkdir -p $output_directory_full/gatk_files	

	wget -O "$output_directory_full/gatk_files/1000G_omni2.5.hg38.vcf.gz" https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/1000G_omni2.5.hg38.vcf.gz
	wget -O "$output_directory_full/gatk_files/1000G_omni2.5.hg38.vcf.gz.tbi" https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/1000G_omni2.5.hg38.vcf.gz.tbi
	wget -O "$output_directory_full/gatk_files/1000G_phase1.snps.high_confidence.hg38.vcf.gz" https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz
	wget -O "$output_directory_full/gatk_files/1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi" https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi
	wget -O "$output_directory_full/gatk_files/Homo_sapiens_assembly38.dbsnp138.vcf" https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf
	bgzip -f $output_directory_full/gatk_files/Homo_sapiens_assembly38.dbsnp138.vcf
	tabix -f $output_directory_full/gatk_files/Homo_sapiens_assembly38.dbsnp138.vcf.gz
	wget -O "$output_directory_full/gatk_files/hapmap_3.3.hg38.vcf.gz" https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/hapmap_3.3.hg38.vcf.gz
	wget -O "$output_directory_full/gatk_files/hapmap_3.3.hg38.vcf.gz.tbi" https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/hapmap_3.3.hg38.vcf.gz.tbi
	

	gatk VariantRecalibrator -R $ref_fa_directory/$ref_fa -V $output_directory_full/${out_prefix}.vcf.gz --resource:hapmap,known=false,training=true,truth=true,prior=15.0 $output_directory_full/gatk_files/hapmap_3.3.hg38.vcf.gz --resource:omni,known=false,training=true,truth=false,prior=12.0 $output_directory_full/gatk_files/1000G_omni2.5.hg38.vcf.gz --resource:1000G,known=false,training=true,truth=false,prior=10.0 $output_directory_full/gatk_files/1000G_phase1.snps.high_confidence.hg38.vcf.gz --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $output_directory_full/gatk_files/Homo_sapiens_assembly38.dbsnp138.vcf.gz -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -mode SNP -O $output_directory_full/${out_prefix}.recal --tranches-file $output_directory_full/${out_prefix}.tranches

	 gatk ApplyVQSR -R $ref_fa_directory/$ref_fa -V $output_directory_full/${out_prefix}.vcf.gz -O $output_directory_full/${out_prefix}.recal.vcf.gz --ts-filter-level 99.0 --tranches-file $output_directory_full/${out_prefix}.tranches --recal-file $output_directory_full/${out_prefix}.recal -mode SNP
 
 
	rm -r $output_directory_full/gatk_files
	

	"""


}



