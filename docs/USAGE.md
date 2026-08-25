# Running PipeVar

This guide explains how PipeVar chooses a workflow and provides examples beyond
the first local run in the [root README](../README.md).

## On this page

- [Choose a route](#how-pipevar-chooses-a-route)
- [Caller selection](#caller-selection)
- [Single-sample examples](#single-sample-examples)
- [Batch examples](#batch-examples)
- [Cluster execution](#cluster-example)
- [Light and targeted analysis](#light-mode)
- [Mitochondrial analysis](#mitochondrial-analysis)
- [Phenotype extraction](#phenotype-extraction)
- [Resume a run](#resume-a-run)

## How PipeVar chooses a route

Four choices determine the main route:

1. input: aligned reads, VCF, or prepared ANNOVAR annotations;
2. sequencing type: `short`, `ont`, or `pacbio` for aligned reads;
3. analysis: small variants, structural/copy-number variants, or combined; and
4. execution: one sample or a sample sheet.

`--mode snp` selects the small-variant branch. `--mode sv` selects the
structural/copy-number and repeat branch. Omitting `--mode` selects the supported
combined route for aligned reads.

| Input | Sequencing | Mode | Main result |
| --- | --- | --- | --- |
| Single BAM/CRAM | `short`, `ont`, `pacbio` | `snp`, `sv`, or omitted | Called and prioritized variants; combined routes add a report |
| Basic sample sheet with BAM/CRAM | `short`, `ont`, `pacbio` | `snp`, `sv`, or omitted | Per-sample batch results |
| Typed sheet with `bam_ngs`/`cram_ngs` | `short`, `ont`, `pacbio` | `snp`, `sv`, or omitted | Per-sample batch results; `--type` selects the route |
| Single VCF | Not applicable | `snp` or `sv` required | Re-annotated and prioritized supplied calls |
| Basic sample sheet with VCF | Not applicable | `snp` or `sv` required | Per-sample VCF analysis |
| Typed sheet with `vcf_snv`/`vcf_sv` | Not applicable | Inferred | Per-sample VCF analysis |
| Prepared ANNOVAR small variants | Not applicable | Small-variant-led | Re-prioritized supplied annotations |
| Typed sheet with `annotated_snv` | Short read when alignment-backed | Small-variant-led or supported batch combination | Prepared evidence with optional structural/copy-number and repeat branches |

Detailed file requirements are in the [input guide](INPUTS.md).

## Caller selection

| Data | Default small-variant caller | Light small-variant caller | Structural/copy-number analysis | Repeat analysis |
| --- | --- | --- | --- | --- |
| Short read | DeepVariant | GATK HaplotypeCaller | Manta and CNVnator | ExpansionHunter |
| Oxford Nanopore | Clair3 | NanoCaller | Sniffles | NanoRepeat |
| PacBio | Clair3 | NanoCaller | Sniffles | NanoRepeat |
| Existing VCF | Supplied calls | Not applicable | Supplied calls | Requires an alignment |

## Single-sample examples

### Short-read combined analysis

Use the complete [first-run command](../README.md#3-run-one-sample-with-docker)
from the root README. The examples below cover routes that differ from it.

### Oxford Nanopore combined analysis

```bash
nextflow run main.nf \
  -profile local_singularity \
  --bam /data/proband.ont.bam \
  --ref_fa /refs/hg38.fa \
  --note /data/proband.note.txt \
  --type ont \
  --out_prefix proband_ont \
  --output_directory /results/proband_ont
```

Use `--type pacbio` for a PacBio aligned-read input.

### Small variants only

```bash
nextflow run main.nf \
  -profile local_docker \
  --bam /data/proband.bam \
  --ref_fa /refs/hg38.fa \
  --hpo /data/proband.hpo.txt \
  --type short \
  --mode snp \
  --out_prefix proband_snp
```

### Structural/copy-number variants and repeats

```bash
nextflow run main.nf \
  -profile local_docker \
  --bam /data/proband.bam \
  --ref_fa /refs/hg38.fa \
  --hpo /data/proband.hpo.txt \
  --type short \
  --mode sv \
  --out_prefix proband_sv
```

### Existing small-variant VCF

```bash
nextflow run main.nf \
  -profile local_docker \
  --vcf /data/proband.vcf \
  --ref_fa /refs/hg38.fa \
  --mode snp \
  --hpo /data/proband.hpo.txt \
  --out_prefix proband_vcf
```

### Prepared ANNOVAR small variants

```bash
nextflow run main.nf \
  -profile local_docker \
  --annotated_snv yes \
  --annovar_txt /data/proband.hg38_multianno.txt \
  --vcf /data/proband.hg38_multianno.vcf \
  --hpo /data/proband.hpo.txt \
  --out_prefix proband_annotated
```

## Batch examples

### Basic sample sheet

```bash
nextflow run main.nf \
  -profile local_docker \
  --input_csv /data/samples.csv \
  --bam true \
  --note no \
  --ref_fa /refs/hg38.fa \
  --type short \
  --output_directory /results/cohort
```

`--note no` means `note_path` contains HPO files. Omit it when `note_path`
contains clinical notes.

### Typed sample sheet

```bash
nextflow run main.nf \
  -profile local_singularity \
  --input_csv /data/samples.csv \
  --ref_fa /refs/hg38.fa \
  --type short \
  --output_directory /results/cohort
```

Do not add `--bam true` or `--vcf true`; each row's `input_kind` identifies the
input. See [Input guide](INPUTS.md#typed-sample-sheet) for schemas and examples.

## Cluster example

The named SLURM profile uses Singularity:

```bash
nextflow run main.nf \
  -profile slurm_singularity \
  --input_csv /shared/cohort/samples.csv \
  --ref_fa /shared/references/hg38.fa \
  --type ont \
  --output_directory /shared/results/cohort
```

Host resource and model paths must be visible at the same absolute locations on
the submit and compute nodes. See [Installation and setup](INSTALLATION.md).

## Light mode

Set `--light yes` with the documented small-variant-only route:

- short reads use GATK HaplotypeCaller instead of DeepVariant;
- long reads use NanoCaller instead of Clair3; and
- PhenoSV uses its light configuration where applicable.

Use the default full configuration for combined analysis.

## Targeted analysis and gene restriction

`--target yes` builds phenotype-ranked regions from the top Phen2Gene results
and restricts supported small-variant calling or annotation. Change the number
of retained genes with `--phen2gene_filter`; the default is `500`.

Use `--gene` to restrict final prioritization to comma-separated symbols or a
one-symbol-per-line file. Typed `annotated_snv` input does not support targeted
region construction.

## Mitochondrial analysis

Enable mitochondrial analysis with `--mito yes` for BAM/CRAM input. It is
supported with `--mode snp` or an omitted mode, but not with `--mode sv`.
Long-read mitochondrial analysis is not supported with `--light yes`.

Short-read mitochondrial execution requires a FASTA index, sequence dictionary,
and BWA sidecars. Long-read execution requires the FASTA index. Both routes
produce normalized calls, annotations, and a prioritized table separate from
the nuclear prioritized VCF.

```bash
nextflow run main.nf \
  -profile local_docker \
  --bam /data/proband.bam \
  --ref_fa /refs/hg38.fa \
  --hpo /data/proband.hpo.txt \
  --type short \
  --mito yes \
  --out_prefix proband
```

## Phenotype extraction

Direct HPO input uses `--hpo`. Clinical-note input uses `--note` and defaults to
PhenoTagger.

For PhenoGPT2, set `--phenotype_extractor phenogpt2 --GPU yes` and configure the
external model paths described in
[Installation and setup](INSTALLATION.md#phenogpt2-and-gpu-execution).

## Resume a run

Add `-resume` when code, inputs, parameters, containers, and mounted resources
are unchanged:

```bash
nextflow run main.nf \
  -resume \
  -profile local_docker \
  --bam /data/proband.bam \
  --ref_fa /refs/hg38.fa \
  --hpo /data/proband.hpo.txt \
  --type short \
  --out_prefix proband
```

Nextflow options such as `-resume` and `-profile` use one hyphen. PipeVar
parameters use two.

Continue with the [parameter reference](PARAMETERS.md), [workflow guide](WORKFLOW.md),
or [output guide](OUTPUTS.md).
