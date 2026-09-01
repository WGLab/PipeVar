# Canonical command recipes

Use these recipes after resolving the route and completing preflight. Replace every angle-bracket placeholder with an absolute path or deliberate value. If values are not available, label the result as a non-runnable template and list every required substitution. Commands run from the PipeVar repository root.

Prefer `local_docker`, `local_singularity`, or `slurm_singularity` explicitly. Nextflow options such as `-profile` and `-resume` use one hyphen; PipeVar parameters use two.

## Single aligned-read routes

Combined nuclear analysis:

```bash
nextflow run main.nf \
  -profile <PROFILE> \
  --bam <ALIGNMENT> \
  --ref_fa <HG38_FASTA> \
  --hpo <HPO_FILE> \
  --type <short|ont|pacbio> \
  --out_prefix <SAMPLE> \
  --output_directory <OUTPUT_DIR>
```

Omitting `--mode` is intentional. For a clinical note, replace `--hpo <HPO_FILE>` with `--note <NOTE_FILE>`; PhenoTagger is the default extractor.

Small variants only: add `--mode snp`. Structural/copy-number variants and repeats: add `--mode sv`. Do not add both and do not write `--mode combined`.

Expected primary result: combined routes normally provide `*.variant_html_report.html`; branch-only routes start with `*.prio.vcf`. Confirm conditional files with [result-review.md](result-review.md).

## Single existing VCF

```bash
nextflow run main.nf \
  -profile <PROFILE> \
  --vcf <VCF> \
  --ref_fa <HG38_FASTA> \
  --mode <snp|sv> \
  --hpo <HPO_FILE> \
  --out_prefix <SAMPLE> \
  --output_directory <OUTPUT_DIR>
```

Omit `--type`. `--ref_fa` is recommended for a single VCF and required with `--target yes`; if deliberately omitted for a non-targeted single-VCF run, explain that choice. Supplied VCFs do not select the aligned-read combined route.

## Single prepared ANNOVAR SNV

```bash
nextflow run main.nf \
  -profile <PROFILE> \
  --annotated_snv yes \
  --annovar_txt <SAMPLE.hg38_multianno.txt> \
  --vcf <SAMPLE.hg38_multianno.vcf> \
  --hpo <HPO_FILE> \
  --out_prefix <SAMPLE> \
  --output_directory <OUTPUT_DIR>
```

Use only a matching TXT/VCF pair. The canonical command intentionally omits `--mode`; `--mode snp` is also accepted. Do not add `--mode sv` or `--target yes`. Consult [route-selection.md](route-selection.md) and [manifests.md](manifests.md) for documented advanced prepared-input combinations.

## Basic batch sheet

Aligned input:

```bash
nextflow run main.nf \
  -profile <PROFILE> \
  --input_csv <SAMPLES.csv> \
  --bam true \
  --note no \
  --ref_fa <HG38_FASTA> \
  --type <short|ont|pacbio> \
  --output_directory <OUTPUT_DIR>
```

`--note no` means the basic sheet's `note_path` column contains HPO files. Omit it when the column contains clinical notes. Add `--mode snp` or `--mode sv` only when a branch-only run is intended.

For a basic VCF sheet, replace `--bam true` and `--type ...` with `--vcf true --mode <snp|sv>`. Keep `--ref_fa`.

## Typed batch sheet

Aligned `bam_ngs`/`cram_ngs` rows:

```bash
nextflow run main.nf \
  -profile <PROFILE> \
  --input_csv <SAMPLES.csv> \
  --ref_fa <HG38_FASTA> \
  --type <short|ont|pacbio> \
  --output_directory <OUTPUT_DIR>
```

Add a supported run-level mode only when requested. Do not add `--bam true`, `--vcf true`, `--hpo`, or `--note`; row fields supply input and phenotype paths.

For typed `vcf_snv` or `vcf_sv` rows, use the same command without `--type` and without `--mode`; `input_kind` selects the branch. Inspect [manifests.md](manifests.md) before using prepared or alignment-backed annotated rows.

## Mitochondrial analysis

Add `--mito yes` to a supported single or batch aligned-read SNP/combined command. Short-read combined example:

```bash
nextflow run main.nf \
  -profile <PROFILE> \
  --bam <ALIGNMENT> \
  --ref_fa <HG38_FASTA> \
  --hpo <HPO_FILE> \
  --type short \
  --mito yes \
  --out_prefix <SAMPLE> \
  --output_directory <OUTPUT_DIR>
```

For ONT or PacBio, change `--type` and retain the normal non-light caller. Do not combine mitochondrial analysis with `--mode sv`; do not combine long-read mitochondrial analysis with `--light yes`. Read [mitochondrial.md](mitochondrial.md) before launch.

## Clinical-note extraction with PhenoGPT2

Only when the user explicitly selects PhenoGPT2 and the environment is configured, replace the phenotype argument and add the required settings:

```text
--note <NOTE_FILE>
--phenotype_extractor phenogpt2
--GPU yes
--phenogpt2_model_host_path <VERSIONED_MODEL_DIR>
```

Negation processing requires its additional model paths. Do not add these flags to HPO-only runs. Read [runtime-and-resume.md](runtime-and-resume.md).

## Resume

For an eligible unchanged run, add `-resume` as a Nextflow option:

```bash
nextflow run main.nf \
  -resume \
  -profile <PROFILE> \
  <THE_SAME_PIPEVAR_PARAMETERS>
```

Do not represent a changed command with placeholder text as runnable. Reproduce the original complete parameter set and apply the decision rules in [runtime-and-resume.md](runtime-and-resume.md).

## Before returning or running a recipe

- Remove every unresolved placeholder or clearly list what the user must supply.
- Quote paths if they contain shell-significant characters; prefer simple absolute paths.
- Do not include clinical-note contents on the command line.
- Explain why mode is present or intentionally absent.
- State the expected primary result without promising conditional artifacts.
- If the user did not authorize execution, return the command and preflight findings only.
