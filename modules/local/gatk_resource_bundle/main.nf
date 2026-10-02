// Download and validate the hg38 VQSR bundle once per workflow execution.
process GATK_RESOURCE_BUNDLE {
    tag 'hg38-v0'
    label 'process_medium'

    container 'broadinstitute/gatk:4.5.0.0'

    output:
    tuple path('gatk_files/1000G_omni2.5.hg38.vcf.gz'),
    path('gatk_files/1000G_omni2.5.hg38.vcf.gz.tbi'),
    path('gatk_files/1000G_phase1.snps.high_confidence.hg38.vcf.gz'),
    path('gatk_files/1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi'),
    path('gatk_files/Homo_sapiens_assembly38.dbsnp138.vcf.gz'),
    path('gatk_files/Homo_sapiens_assembly38.dbsnp138.vcf.gz.tbi'),
    path('gatk_files/hapmap_3.3.hg38.vcf.gz'),
    path('gatk_files/hapmap_3.3.hg38.vcf.gz.tbi'),
    emit: resources

    script:
    """
    mkdir -p gatk_files

    wget -q -O gatk_files/1000G_omni2.5.hg38.vcf.gz \\
        https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/1000G_omni2.5.hg38.vcf.gz
    wget -q -O gatk_files/1000G_omni2.5.hg38.vcf.gz.tbi \\
        https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/1000G_omni2.5.hg38.vcf.gz.tbi
    wget -q -O gatk_files/1000G_phase1.snps.high_confidence.hg38.vcf.gz \\
        https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz
    wget -q -O gatk_files/1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi \\
        https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi
    wget -q -O gatk_files/Homo_sapiens_assembly38.dbsnp138.vcf \\
        https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf
    wget -q -O gatk_files/hapmap_3.3.hg38.vcf.gz \\
        https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/hapmap_3.3.hg38.vcf.gz
    wget -q -O gatk_files/hapmap_3.3.hg38.vcf.gz.tbi \\
        https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/hapmap_3.3.hg38.vcf.gz.tbi

    cat <<-'END_CHECKSUMS' | md5sum --check --strict
    345b811a791f6704996b8bca4c8a16d9  gatk_files/1000G_omni2.5.hg38.vcf.gz
    37150862f34b237763e737f840fb24e5  gatk_files/1000G_omni2.5.hg38.vcf.gz.tbi
    b2979b47800b59b41920bf5432c4b2a0  gatk_files/1000G_phase1.snps.high_confidence.hg38.vcf.gz
    1610d5349edae9df29bd0ce62fdea8ce  gatk_files/1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi
    d05ac6b9a247a21ce0030c7494194da9  gatk_files/hapmap_3.3.hg38.vcf.gz
    00897695134755db73b7a360244a4cf1  gatk_files/hapmap_3.3.hg38.vcf.gz.tbi
    END_CHECKSUMS

    test "\$(stat -c '%s' gatk_files/Homo_sapiens_assembly38.dbsnp138.vcf)" = '10950827213'
    bgzip -f gatk_files/Homo_sapiens_assembly38.dbsnp138.vcf
    tabix -f -p vcf gatk_files/Homo_sapiens_assembly38.dbsnp138.vcf.gz
    """

    stub:
    """
    mkdir -p gatk_files
    for file in \\
        1000G_omni2.5.hg38.vcf.gz \\
        1000G_phase1.snps.high_confidence.hg38.vcf.gz \\
        Homo_sapiens_assembly38.dbsnp138.vcf.gz \\
        hapmap_3.3.hg38.vcf.gz
    do
        echo '' | gzip > "gatk_files/\$file"
        touch "gatk_files/\$file.tbi"
    done
    """
}
