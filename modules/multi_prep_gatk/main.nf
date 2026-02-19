
// Prepare batch GATK resource files and bundle inputs for HaplotypeCaller.
process multi_prep_gatk {
	container ='broadinstitute/gatk:4.5.0.0'

        input:
	tuple val(out_prefix), path(bam), path(index_file)

	output:
	tuple val(out_prefix), path(bam), path(index_file), path("gatk_files/1000G_omni2.5.hg38.vcf.gz"), path("gatk_files/1000G_omni2.5.hg38.vcf.gz.tbi"), path("gatk_files/1000G_phase1.snps.high_confidence.hg38.vcf.gz"), path("gatk_files/1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi"), path("gatk_files/Homo_sapiens_assembly38.dbsnp138.vcf.gz"), path("gatk_files/Homo_sapiens_assembly38.dbsnp138.vcf.gz.tbi"), path("gatk_files/hapmap_3.3.hg38.vcf.gz"), path("gatk_files/hapmap_3.3.hg38.vcf.gz.tbi") 

	script:

	"""


	
	mkdir -p gatk_files	

	wget -O "./gatk_files/1000G_omni2.5.hg38.vcf.gz" https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/1000G_omni2.5.hg38.vcf.gz
	wget -O "./gatk_files/1000G_omni2.5.hg38.vcf.gz.tbi" https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/1000G_omni2.5.hg38.vcf.gz.tbi
	wget -O "./gatk_files/1000G_phase1.snps.high_confidence.hg38.vcf.gz" https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz
	wget -O "./gatk_files/1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi" https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi
	wget -O "./gatk_files/Homo_sapiens_assembly38.dbsnp138.vcf" https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf
	bgzip -f ./gatk_files/Homo_sapiens_assembly38.dbsnp138.vcf
	tabix -f ./gatk_files/Homo_sapiens_assembly38.dbsnp138.vcf.gz
	wget -O "./gatk_files/hapmap_3.3.hg38.vcf.gz" https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/hapmap_3.3.hg38.vcf.gz
	wget -O "./gatk_files/hapmap_3.3.hg38.vcf.gz.tbi" https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/hapmap_3.3.hg38.vcf.gz.tbi
	

 
 
	

	"""


}



