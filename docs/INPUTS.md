# PipeVar input guide

PipeVar accepts aligned reads, existing variant calls, prepared ANNOVAR
annotations, and phenotype information. This page describes the required files
for single-sample and sample-sheet execution.

## On this page

- [Input concepts](#input-concepts)
- [References and indexes](#reference-files)
- [Phenotype files](#phenotype-files)
- [Single-file inputs](#single-aligned-read-input)
- [Choose a sample sheet](#choose-a-sample-sheet)
- [Basic sample sheet](#basic-sample-sheet)
- [Typed sample sheet](#typed-sample-sheet)
- [Generate a sample sheet](#generate-a-sample-sheet)

## Input concepts

| Term | Meaning |
| --- | --- |
| BAM or CRAM | Reads already aligned to a reference genome |
| VCF | Existing variant calls |
| ANNOVAR multianno | Prepared hg38 annotation TXT and matching VCF |
| HPO | Human Phenotype Ontology terms describing the patient's findings |
| Clinical note | Free text converted to HPO terms by a configured extractor |
| Sample sheet | CSV file describing multiple inputs, one row per sample |

## Reference files

Alignment-driven workflows require `--ref_fa <FASTA>`. Any supplied FASTA must
have a samtools index at `<reference>.fai`.

PipeVar's active ANNOVAR modules operate on hg38. Input coordinates, aligned
reads, FASTA sequence, indexes, and annotation databases must use the same
reference build.

Short-read mitochondrial analysis additionally requires:

- `<reference-basename>.dict` beside the FASTA; and
- BWA sidecars `<reference>.amb`, `.ann`, `.bwt`, `.pac`, and `.sa`.

Long-read mitochondrial analysis requires the FASTA index but not the BWA
sidecars or sequence dictionary.

## Alignment indexes

- BAM input requires `<file>.bam.bai` or, where discovered by the sample-sheet
  helper, `<file-basename>.bai`.
- CRAM input requires `<file>.cram.crai` or `<file-basename>.crai`.
- CRAM decoding requires the exact matching reference FASTA.

## Phenotype files

Every analyzed sample requires one phenotype source:

- `--hpo <FILE>` for HPO terms; or
- `--note <FILE>` for a clinical note.

Clinical notes use PhenoTagger by default. PhenoGPT2 configuration is described
in [Installation and setup](INSTALLATION.md#phenogpt2-and-gpu-execution).

HPO files should contain one HPO identifier per line. Keep identifiers in their
standard `HP:NNNNNNN` form.

## Single aligned-read input

Required arguments are:

```text
--bam <BAM_OR_CRAM>
--ref_fa <HG38_FASTA>
(--hpo <HPO_FILE> | --note <CLINICAL_NOTE>)
--type <short|ont|pacbio>
```

Use `--mode snp` for small variants, `--mode sv` for structural/copy-number
analysis, or omit `--mode` for the supported combined analysis.

## Single VCF input

Use an explicit mode:

```text
--vcf <VCF>
--mode <snp|sv>
(--hpo <HPO_FILE> | --note <CLINICAL_NOTE>)
```

Supplying `--ref_fa` is recommended and is required for targeted region
construction. VCF sample-sheet routes require it.

## Prepared ANNOVAR input

For a single prepared small-variant input, provide the paired multianno files:

```text
--annotated_snv yes
--annovar_txt sample.hg38_multianno.txt
--vcf sample.hg38_multianno.vcf
(--hpo <HPO_FILE> | --note <CLINICAL_NOTE>)
```

PipeVar validates and passes the original pair downstream without publishing
full-size validated copies.

For advanced batch combinations with an annotated structural-variant VCF or an
alignment, use a typed sample sheet.

## Choose a sample sheet

| Format | Best for | How input is selected |
| --- | --- | --- |
| Basic sample sheet | One input family shared by every row | Launch with `--bam true` or `--vcf true` |
| Typed sample sheet | Explicit or advanced per-row paths | Each row supplies `input_kind` and its matching path fields |

## Basic sample sheet

The basic sample sheet uses one path column for the selected input family:

```csv
sample,file_path,note_path
sample1,/data/sample1.bam,/phenotypes/sample1.hpo.txt
```

The column is named `note_path`, but it can contain either a clinical note or an
HPO file. This example uses HPO input and therefore launches with `--note no`.

| Column | Required | Meaning |
| --- | --- | --- |
| `sample` | Yes | Unique sample identifier and output prefix |
| `file_path` | Yes | BAM, CRAM, or VCF path selected by the launch command |
| `note_path` | Yes | Clinical note by default, or HPO file with `--note no` |
| `age_of_onset` | No | Integer or integer plus `d`, `m`, or `y`; bare integers become years |
| `age` | No | Alternate age column used only when `age_of_onset` is absent |

Launch aligned reads with:

```bash
nextflow run main.nf \
  -profile local_docker \
  --input_csv /data/samples.csv \
  --bam true \
  --note no \
  --ref_fa /refs/hg38.fa \
  --type short
```

Launch VCF input with `--vcf true --mode <snp|sv>` instead of `--bam true`.

## Typed sample sheet

The typed sheet describes the input in each row. `input_kind` selects the input
family, and the matching path column points to the file.

A common aligned-read example is:

```csv
sample,input_kind,phenotype_path,phenotype_format,age_of_onset,snv_txt_path,snv_vcf_path,sv_vcf_path,vcf_path,alignment_path,alignment_index_path
sample1,bam_ngs,/phenotypes/sample1.hpo.txt,hpo,,,,,,/data/sample1.bam,/data/sample1.bam.bai
```

For this row, the meaningful values are `sample`, `input_kind`, the phenotype
path/format, and the alignment path/index. Empty fields retain the exact shared
header but are unused for this input kind.

Required base columns are:

| Column | Values | Meaning |
| --- | --- | --- |
| `sample` | Unique text | Join key and output prefix |
| `input_kind` | See below | Input family for the row |
| `phenotype_path` | File path | Clinical note or HPO file |
| `phenotype_format` | `clinical_note` or `hpo` | How to interpret the phenotype file |

Conditional path columns are:

| `input_kind` | Main path columns |
| --- | --- |
| `bam_ngs`, `cram_ngs` | `alignment_path`; optional explicit `alignment_index_path` |
| `vcf_snv`, `vcf_sv` | `vcf_path` |
| `annotated_snv` | `snv_txt_path` and `snv_vcf_path`; optional `sv_vcf_path`, `alignment_path`, and `alignment_index_path` |

`age_of_onset` is optional. The remaining path columns may be present with blank
values when they do not apply to a row.

Typed-sheet constraints are:

- sample identifiers must be unique;
- rows must use a supported homogeneous input-kind group, except BAM and CRAM
  aligned-read rows may be combined;
- all `annotated_snv` rows must consistently provide or omit an alignment;
- all `annotated_snv` rows must consistently provide or omit `sv_vcf_path`;
- non-annotated sheets use one phenotype format per run; and
- alignment-backed `annotated_snv` rows currently support short reads.

Prepared-annotation batch layouts support:

- annotated small variants only;
- annotated small and structural variants without an alignment;
- annotated small variants plus an alignment for structural-variant and repeat
  analysis; or
- annotated small and structural variants plus an alignment.

Routes without an alignment cannot perform read-backed repeat or mitochondrial
analysis.

Launch the typed sheet without `--bam true` or `--vcf true`:

```bash
nextflow run main.nf \
  -profile local_docker \
  --input_csv /data/samples.csv \
  --ref_fa /refs/hg38.fa \
  --type short
```

For `vcf_snv` and `vcf_sv`, the mode is inferred from `input_kind`; omit an
explicit contradictory `--mode`.

## Generate a sample sheet

Run the interactive helper:

```bash
./scripts/generate_input_csv.sh
```

It can create either sample-sheet format and add `sv_vcf_path` to an existing
typed annotated-SNV sheet.

Review generated files before launching because:

- only the selected directory's top level is scanned;
- matching uses exact, normalized, containment, and first-token strategies;
- ambiguous or incomplete pairs can be skipped;
- a missing index can produce a blank index field;
- commas and quoted CSV fields are not supported; and
- duplicate inferred sample prefixes keep the first match with a warning.

The update action requires a separate output path and does not overwrite its
source sheet.

Continue with [Running PipeVar](USAGE.md) or return to the
[documentation map](../README.md#documentation).
