FROM python:3.11-slim

ENV MITO_DB_ROOT=/opt/mito_db \
    MITO_DB_RAW=/opt/mito_db/raw \
    MITO_DB_PREPARED=/opt/mito_db/prepared \
    PATH=/opt/mito/bin:$PATH

WORKDIR /opt/mito

COPY requirements.txt /opt/mito/requirements.txt
RUN pip install --no-cache-dir -r /opt/mito/requirements.txt

COPY mitomap_variant_evidence.tsv /opt/mito_db/raw/mitomap_variant_evidence.tsv
COPY mitotip_scores.txt /opt/mito_db/raw/mitotip_scores.txt
COPY t-APOGEE_2024.0.1.txt /opt/mito_db/raw/t-APOGEE_2024.0.1.txt
COPY MitImpact_db_3.1.3.txt /opt/mito_db/raw/MitImpact_db_3.1.3.txt
COPY build_mito_manifest.py /opt/mito/bin/build_mito_manifest.py
COPY prepare_mito_db.py /opt/mito/bin/prepare_mito_db.py
COPY annotate_mito_variants.py /opt/mito/bin/annotate_mito_variants.py
COPY prioritize_mito_variants.py /opt/mito/bin/prioritize_mito_variants.py

RUN chmod +x /opt/mito/bin/*.py \
    && mkdir -p /opt/mito_db/prepared \
    && python3 /opt/mito/bin/prepare_mito_db.py --raw-dir /opt/mito_db/raw --out-dir /opt/mito_db/prepared \
    && python3 /opt/mito/bin/build_mito_manifest.py --raw-dir /opt/mito_db/raw --prepared-dir /opt/mito_db/prepared --out /opt/mito_db/manifest.tsv

WORKDIR /work
