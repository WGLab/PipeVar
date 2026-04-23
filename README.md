# PipeVar_annotated_snv

PipeVar_annotated_snv is a Nextflow DSL2 workflow for rare-disease variant prioritization from short-read and long-read data, with support for a standalone pre-annotated ANNOVAR SNP input mode.
It preserves the existing BAM/CRAM, VCF, SV, STR, mito, and SCRAMBLE analysis paths from `PipeVar_mito` and adds an ANNOVAR TXT + VCF SNP-prioritization path that skips ANNOVAR runtime re-annotation.

## What PipeVar does

- Calls and prioritizes SNP/indel variants.
- Accepts pre-annotated ANNOVAR SNP input (`.hg38_multianno.txt` + `.hg38_multianno.vcf`) and prioritizes it without rerunning ANNOVAR.
- Calls and prioritizes structural variants (SV).
- Runs repeat expansion analysis (short-read and long-read paths).
- Optionally calls and annotates mitochondrial variants with Mutect2 plus bundled mtDNA evidence databases.
- Uses phenotype inputs (`--hpo` or clinical note via `--note`) for phenotype-guided ranking.
- Supports single-sample mode and CSV batch mode.

## Runtime model

PipeVar is designed for containerized execution.

- Supported container backends:
  - Singularity
  - Docker
- Tested/primary scheduler profile:
  - SLURM (`standard` / `slurm_singularity`)
- Also available:
  - local executor with Singularity
  - local executor with Docker

## Execution profiles

Defined in `nextflow.config`:

- `standard`
  - SLURM + Singularity (default profile behavior)
- `slurm_singularity`
  - Explicit SLURM + Singularity
- `local_singularity`
  - Local executor + Singularity
- `local_docker`
  - Local executor + Docker

All Singularity/Docker profiles mount:

- `--annovar_host_path` -> `/annovar`
- `--phenosv_host_path` -> `/PhenoSV/train_data`

Mitochondrial databases are baked into the mito annotation Docker image and do not require extra runtime bind mounts.

## Setup

### 1) Clone repository

```bash
git clone https://github.com/WGLab/PipeVar.git
cd PipeVar
```

### 2) External data/software prerequisites

PipeVar_annotated_snv expects ANNOVAR and PhenoSV resources to be available for legacy re-annotation modes.
The new `--annotated_snv yes` mode does not require ANNOVAR runtime databases because it consumes pre-annotated ANNOVAR outputs directly.

ANNOVAR registration/download:

- https://www.openbioinformatics.org/annovar/annovar_download_form.php

Then run setup script:

```bash
# Full setup
./setup.sh

# Light PhenoSV setup
./setup.sh light
```

By default, setup expects:

- ANNOVAR at `./annovar`
- PhenoSV resources downloaded under `./PhenoSV_model`

You can override both locations (recommended for HPC/shared filesystems):

```bash
./setup.sh --annovar-dir=/shared/apps/annovar --phenosv-dir=/shared/data/PhenoSV_model
```

The setup script prepares required assets and writes host-path references used by runtime mounts.
It now also writes a local override file, `.pipevar.user.config`, with:

- a persisted default execution profile (`manifest.defaultProfile`)
- persisted bind source paths:
  - `params.annovar_host_path`
  - `params.phenosv_host_path`

So after setup, users can run without repeatedly passing `-profile` and bind-path params.

Non-interactive setup example:

```bash
./setup.sh --non-interactive --profile=local_docker \\
  --annovar-dir=/data/annovar \\
  --phenosv-dir=/data/PhenoSV_model \\
  --annovar-bind=/data/annovar \\
  --phenosv-bind=/data/PhenoSV_model
```

## Pre-annotated ANNOVAR SNP mode

PipeVar_annotated_snv adds a standalone SNP-only input mode for pre-annotated ANNOVAR output:

- enable with `--annotated_snv yes`
- required inputs:
  - `--annovar_txt sample.hg38_multianno.txt`
  - `--vcf sample.hg38_multianno.vcf`
  - one phenotype source: `--note <FILE>` or `--hpo <FILE>`
- this mode skips the `ANNOVAR` module and reuses the existing `Phen2gene`, `RankScore`, `RankVar`, and `snp_prio` stack
- if you also provide BAM/CRAM input (`--bam ... --type short --ref_fa ...`), the same annotated-SNV mode now runs short-read SV, CNV, STR, and optional mito analysis alongside the imported SNP evidence
- `--mode sv` is rejected
- `--target yes` is rejected
- `--ref_fa` is not required for SNP-only annotated-SNV runs
- `--ref_fa` is required for annotated-SNV + BAM/CRAM hybrid runs

Example:

```bash
nextflow run main.nf \
  -profile local_docker \
  --annotated_snv yes \
  --annovar_txt /data/sample.hg38_multianno.txt \
  --vcf /data/sample.hg38_multianno.vcf \
  --hpo /data/sample.hpo.txt \
  --out_prefix sample1
```

Hybrid short-read example:

```bash
nextflow run main.nf \
  -profile local_docker \
  --annotated_snv yes \
  --annovar_txt /data/sample.hg38_multianno.txt \
  --vcf /data/sample.hg38_multianno.vcf \
  --bam /data/sample.cram \
  --ref_fa /refs/hg38.fa \
  --type short \
  --mito yes \
  --hpo /data/sample.hpo.txt \
  --out_prefix sample1_all_ngs
```

## Mitochondrial analysis

PipeVar_annotated_snv keeps the opt-in mitochondrial branch for short-read BAM/CRAM input:

- enable with `--mito yes`
- supported only with `--type short`
- supported only with BAM/CRAM input, not VCF-only mode
- supported with `--mode snp` or when `--mode` is omitted

The mito branch emits separate outputs and does not modify the existing nuclear `.prio.vcf` outputs.

Bundled mito database sources in the annotation image:

- `MITOMAP`
- `MitoTip`
- `t-APOGEE`
- `MitImpact`

Example:

```bash
nextflow run main.nf \
  -profile local_docker \
  --bam sample.bam \
  --ref_fa ref.fa \
  --hpo sample.hpo.txt \
  --type short \
  --mode snp \
  --mito yes \
  --out_prefix sample1
```

## SCRAMBLE mobile-element analysis

PipeVar_annotated_snv also supports an opt-in SCRAMBLE branch for nuclear short-read SV/MEI analysis:

- enable with `--scramble yes`
- supported only with `--type short`
- supported only with BAM/CRAM input, not VCF-only mode
- supported with `--mode sv` or when `--mode` is omitted
- executed per sample in both single-sample and CSV batch BAM/CRAM modes

PipeVar runs SCRAMBLE as one internal short-read subworkflow step that performs:

- soft-clipped cluster identification
- cluster analysis into MEI/deletion outputs

The shared SCRAMBLE image is expected to provide:

- `cluster_identifier`
- `cluster_analysis/bin/SCRAMble.R`
- `cluster_analysis/resources/MEI_consensus_seqs.fa`

Important notes:

- SCRAMBLE is intended for short-read WGS MEI discovery/genotyping.
- BAM and CRAM must be indexed.
- CRAM runs require the matching reference FASTA.
- PipeVar now builds the reference BLAST database SCRAMBLE needs from `--ref_fa`, so users do not need to prepare `.nhr/.nin/.nsq` sidecar files ahead of time.
- SCRAMBLE CRAM support remains provisional because upstream `cluster_identifier` reference resolution still needs dedicated validation.
- SCRAMBLE download/build is handled outside this repo; the shared image should already contain the bundled assets.

## Input modes

## Single-sample BAM/CRAM mode

Required:

- `--bam <FILE>`
- `--ref_fa <FILE>`
- one phenotype source:
  - `--note <FILE>` (clinical note; PipeVar runs PhenoTagger)
  - `--hpo <FILE>` (HPO term file)

Optional:

- `--mode <snp|sv>` to run only one branch

## Single-sample VCF mode

Required:

- `--vcf <FILE>`
- `--ref_fa <FILE>`
- `--mode <snp|sv>`
- one phenotype source (`--note` or `--hpo`)

## CSV batch mode (Unified manifest)

`PipeVar_annotated_snv` now supports a unified batch manifest.

Required common columns:

- `sample,input_kind,phenotype_path,phenotype_format`

Conditional file columns:

- `snv_txt_path,snv_vcf_path,vcf_path,alignment_path,alignment_index_path`

Supported `input_kind` values:

- `annotated_snv`
- `vcf_snv`
- `vcf_sv`
- `bam_ngs`
- `cram_ngs`

Supported `phenotype_format` values:

- `clinical_note`
- `hpo`

Example annotated-SNV row:

```csv
sample,input_kind,phenotype_path,phenotype_format,age_of_onset,snv_txt_path,snv_vcf_path,vcf_path,alignment_path,alignment_index_path
P001,annotated_snv,/data/pheno/P001.hpo.txt,hpo,7y,/data/annovar/P001.hg38_multianno.txt,/data/annovar/P001.hg38_multianno.vcf,,,
```

Hybrid annotated-SNV + alignment row:

```csv
sample,input_kind,phenotype_path,phenotype_format,age_of_onset,snv_txt_path,snv_vcf_path,vcf_path,alignment_path,alignment_index_path
P001,annotated_snv,/data/pheno/P001.hpo.txt,hpo,7y,/data/annovar/P001.hg38_multianno.txt,/data/annovar/P001.hg38_multianno.vcf,,/data/ngs/P001.cram,/data/ngs/P001.cram.crai
```

The helper script `scripts/generate_input_csv.sh` now writes this unified schema.

## CSV batch mode (BAM/CRAM, legacy)

Required:

- `--input_csv <FILE>`
- `--bam true`
- `--ref_fa <FILE>`

Expected CSV columns:

- `sample,file_path,note_path`
- Optional age column for CSV prioritization flows:
  - `sample,file_path,note_path,age_of_onset`
  - `sample,file_path,note_path,age`
  - If both are present, `age_of_onset` is used.
  - Age is interpreted per row (per sample), not globally.
  - Empty age is allowed and treated as not provided.
  - Non-empty age must be `xd`/`xm`/`xy` or integer years.
  - Examples: `10d`, `9m`, `7y`, `7` (`7` is normalized to `7y`).

Phenotype handling in CSV mode:

- default: `note_path` is treated as clinical note (PhenoTagger ON)
- if `--note no`: `note_path` is treated as HPO file (PhenoTagger OFF)

## CSV batch mode (VCF, legacy)

Required:

- `--input_csv <FILE>`
- `--vcf true`
- `--ref_fa <FILE>`
- `--mode <snp|sv>`

Expected CSV columns:

- `sample,file_path,note_path`
- Optional age column for CSV prioritization flows:
  - `sample,file_path,note_path,age_of_onset`
  - `sample,file_path,note_path,age`
  - If both are present, `age_of_onset` is used.
  - Age is interpreted per row (per sample), not globally.
  - Empty age is allowed and treated as not provided.
  - Non-empty age must be `xd`/`xm`/`xy` or integer years.
  - Examples: `10d`, `9m`, `7y`, `7` (`7` is normalized to `7y`).

## Core parameters

- `--bam <FILE>`: single BAM/CRAM input (mutually exclusive with `--vcf` in single-file mode)
- `--vcf <FILE>`: single VCF input
- `--input_csv <FILE>`: manifest for batch processing
- `--ref_fa <FILE>`: reference FASTA
- `--out_prefix <STRING>`: output prefix (single-sample mode)
- `--output_directory <DIR>`: publish directory (default: launch directory)
- `--mode <snp|sv>`: restrict to SNP or SV branch
- `--type <ont|pacbio|short>`: sequencing type for BAM/CRAM flows
- `--light <yes|no>`: enable lightweight models/callers where supported
- `--genome <hg38|grch38>`: genome build for ExpansionHunter catalog selection
- `--target <yes|no>`: restrict SNP calling to phenotype-derived gene BED
- `--scramble <yes|no>`: enable SCRAMBLE MEI calling in short-read SV/all-NGS BAM/CRAM paths (default: `no`)
- `--phen2gene_filter <INT>`: top-N genes retained for targeted mode (default: 500)
- `--rankscore <FLOAT>`: RankScore threshold (default: 0.50)
- `--rankscore_softwares <CSV>`: comma-separated RankScore software names for score aggregation (default: all built-in tools)
- `--gnomad <FLOAT>`: max AF threshold for SNP prioritization (default: 0.0001)
- `--inheritance_mode <ml|omim|gnomad>`: inheritance assignment backend for prioritization (default: `ml`)
- `--include_clinvar_report <yes|no>`: include ClinVar-only calls in final prioritized reports (default: `yes`)
- `--allow_unphased_comphet <yes|no>`: allow unphased `0/1` or `1/0` AR pairs as compound het in final prioritization (default: `no`)
- `--gq <INT>`: genotype quality threshold (default: 20)
- `--ad <INT>`: allele depth threshold (default: 15)
- `--note <FILE|no>`: phenotype note input, or `no` in CSV mode to interpret `note_path` as HPO file
- `--hpo <FILE>`: phenotype HPO file
- `--help`: print help

## Important behavior updates

### Unified light behavior for SNP/all workflows

`--light yes` no longer requires separate SNP/all workflow selection in `main.nf`.
The workflow now uses unified subworkflows and switches SNP caller internally by mode:

- short-read SNP caller:
  - default: `deepvariant`
  - `--light yes`: `haplotypecaller`
- long-read SNP caller:
  - default: `clair3`
  - `--light yes`: `nanocaller`

`--light yes` also enables PhenoSV-light model through config (`ext.args`).

### ExpansionHunter catalog selection

Catalog path is selected from `--genome` for both single and batch modes:

- `hg38` -> `/hg38/variant_catalog.json`
- `grch38` -> `/EH_grch38/variant_catalog.json`

## Example commands

### Single-sample long-read full analysis

```bash
nextflow run main.nf \
  -profile standard \
  --bam /data/p1.bam \
  --ref_fa /refs/hg38.fa \
  --note /data/p1_note.txt \
  --out_prefix p1 \
  --type ont
```

### Single-sample short-read full analysis (light)

```bash
nextflow run main.nf \
  -profile standard \
  --bam /data/p2.bam \
  --ref_fa /refs/hg38.fa \
  --hpo /data/p2_hpo.txt \
  --out_prefix p2 \
  --type short \
  --light yes
```

### Single-sample short-read SV analysis with SCRAMBLE

```bash
nextflow run main.nf \
  -profile standard \
  --bam /data/p2.bam \
  --ref_fa /refs/hg38.fa \
  --hpo /data/p2_hpo.txt \
  --out_prefix p2_sv \
  --type short \
  --mode sv \
  --scramble yes
```

### Single-sample VCF SNP re-annotation/prioritization

```bash
nextflow run main.nf \
  -profile local_docker \
  --vcf /data/p3.vcf \
  --mode snp \
  --ref_fa /refs/hg38.fa \
  --hpo /data/p3_hpo.txt \
  --rankscore_softwares "REVEL,AlphaMissense,CADD_raw" \
  --out_prefix p3
```

### CSV batch BAM mode with HPO file in `note_path`

```bash
nextflow run main.nf \
  -profile slurm_singularity \
  --input_csv /data/samples.csv \
  --bam true \
  --note no \
  --ref_fa /refs/hg38.fa \
  --type short
```

### CSV batch short-read SV mode with SCRAMBLE

```bash
nextflow run main.nf \
  -profile slurm_singularity \
  --input_csv /data/samples.csv \
  --bam true \
  --note no \
  --ref_fa /refs/hg38.fa \
  --type short \
  --mode sv \
  --scramble yes
```

### CSV batch VCF mode (SV only)

```bash
nextflow run main.nf \
  -profile local_singularity \
  --input_csv /data/sv_samples.csv \
  --vcf true \
  --mode sv \
  --ref_fa /refs/hg38.fa
```

## Expected outputs (high-level)

Outputs are published to `--output_directory`.
Exact files depend on `--mode`, `--type`, and input type.

### SNP-related outputs

- caller output (depends on type/light):
  - `*.deepvariant.vcf.gz` (short default)
  - `*.recal.vcf.gz` (short light / HaplotypeCaller path)
  - `*.clair3.vcf.gz` (long default)
  - `*.nanocaller.vcf.gz` (long light)
- annotation/prioritization:
  - `*.clinvar.txt`
  - `*.rank_var.tsv`
  - `*.rankscore_filtered.tsv`
  - ANNOVAR intermediate/final files (`*.hg38_multianno.*`)

### SV-related outputs

- short-read SV:
  - `*_manta.vcf`
  - `*_scramble.vcf` when `--scramble yes`
  - `*.shortread_sv.merged.vcf` when multiple short-read SV/MEI callers are merged
- long-read SV:
  - `*.sniffles.vcf.gz`
- downstream SV prioritization:
  - `*.exonic.vcf`
  - `*.phenosv.filtered.tsv` (or corresponding filtered artifacts)

### Repeat expansion outputs

- short-read:
  - `*.json` (ExpansionHunter raw output)
  - `*.eh.tsv` (filtered disease-threshold loci)
- long-read:
  - NanoRepeat result files (`*_nanoRepeat_output.tsv`, related summary files)

### Phenotype intermediate outputs

- `*_phenotagger_patient_hpo.txt`
- Phen2Gene ranking outputs (`*_phen2gene*`)

## Resource/retry behavior

Configured in `nextflow.config`:

- global process retry strategy:
  - `errorStrategy = 'retry'`
  - `maxRetries = 3`
- CPU/memory/time vary by process via `withName` blocks.

## Notes and pitfalls

- `--input_csv` requires either `--bam true` or `--vcf true`.
- In single-file mode, at least one of `--note <FILE>` or `--hpo <FILE>` is required.
- For single VCF mode, `--mode` must be provided.
- Reference index (`.fai`) must exist.
- BAM/CRAM index must exist (`.bai`/`.crai`) for alignment-driven paths.
- If using Singularity/Docker profiles, ensure `--annovar_host_path` and `--phenosv_host_path` point to valid host locations.

## Software/components used

### SNP calling

- DeepVariant
- GATK HaplotypeCaller (+ VQSR flow in relevant path)
- Clair3
- NanoCaller

### SV calling/prioritization

- Sniffles
- Manta
- SURVIVOR
- PhenoSV
- ANNOVAR SV annotation module

### Repeat expansion

- ExpansionHunter
- NanoRepeat

### Annotation/ranking/phenotype

- ANNOVAR
- RankVar
- RankScore filtering path
- Phen2Gene
- PhenoTagger
- Longphase prioritization helpers

## Status

PipeVar is under active development. If behavior seems inconsistent with this README,
`main.nf` help output and `nextflow.config` are the source of truth.
