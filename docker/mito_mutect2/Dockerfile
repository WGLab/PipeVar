FROM broadinstitute/gatk:4.5.0.0

USER root
RUN apt-get update \
    && apt-get install -y --no-install-recommends bwa samtools bcftools tabix \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /work
