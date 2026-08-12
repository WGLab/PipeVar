# PipeVar_mito

PipeVar_mito is a Nextflow DSL2 workflow for rare-disease variant prioritization from short-read and long-read data. It keeps the nuclear SNV/indel, SV/CNV/MEI, repeat-expansion, phenotype-ranking, and reporting paths from PipeVar, and adds an opt-in mitochondrial analysis branch for BAM/CRAM inputs.

## Overview

PipeVar_mito supports:

- SNV/indel calling, annotation, and prioritization.
- Structural variant and CNV calling, annotation, and prioritization.
- Short-read and long-read repeat-expansion analysis.
- Optional short-read mobile-element analysis with xTEA.
- Optional mitochondrial analysis with Mutect2 for short reads or Clair3 plus postprocessing for long reads.
- Phenotype-guided ranking from HPO files or clinical notes processed through PhenoTagger.
- Single-sample and CSV batch execution.

## Quick Start

Clone the repository:

```bash
git clone https://github.com/WGLab/PipeVar.git PipeVar_mito
cd ./PipeVar_mito
```

Run setup after preparing ANNOVAR and PhenoSV resources:

```bash
./setup.sh --annovar-dir=/shared/apps/annovar --phenosv-dir=/shared/data/PhenoSV_model
```

Minimal single-sample BAM/CRAM run:

```bash
nextflow run main.nf \
  -profile standard \
  --bam /data/sample.bam \
  --ref_fa /refs/hg38.fa \
  --hpo /data/sample.hpo.txt \
  --type ont \
  --out_prefix sample1
```

Minimal legacy CSV batch run with HPO files in `note_path`:

```bash
nextflow run main.nf \
  -profile slurm_singularity \
  --input_csv /data/samples.csv \
  --bam true \
  --note no \
  --ref_fa /refs/hg38.fa \
  --type short
```

## Installation And Runtime

PipeVar_mito is designed for containerized execution. Execution profiles are defined in `nextflow.config`.

| Profile | Executor | Container backend |
| --- | --- | --- |
| `standard` | SLURM | Singularity |
| `slurm_singularity` | SLURM | Singularity |
| `local_singularity` | local | Singularity |
| `local_docker` | local | Docker |

All Singularity/Docker profiles mount these host paths:

| Parameter | Container path | Purpose |
| --- | --- | --- |
| `--annovar_host_path` | `/annovar` | ANNOVAR installation and databases |
| `--phenosv_host_path` | `/PhenoSV/train_data` | PhenoSV model resources |

Mitochondrial annotation databases are baked into the mito annotation image and do not require extra runtime bind mounts.

PhenoGPT2 weights are deliberately not included in its image. When a clinical
note is routed to PhenoGPT2, provide a complete, immutable, versioned checkpoint
directory with `--phenogpt2_model_host_path`. That directory is mounted
read-only only for the PhenoGPT2 task. The default generated-model cache is
task-local; an optional pre-created `--phenogpt2_cache_host_path` is mounted
read-write only for that task. DeepVariant never receives either mount.
The `phenogpt2_0.3-depslim` commands below apply after that candidate has passed
the documented GPU gates and its immutable digest has been pinned in both
PhenoGPT2 modules.

Release status (2026-07-23): the pinned local builds measured 23.36 GB for
`artifact-slim`, 20.65 GB for `dependency-slim`, and 14.69 GB for
the compiler-free `runtime-slim` candidate (Docker uncompressed image size).
Cold-cache vLLM initialization showed that Triton requires a native compiler at
runtime. The corrected runtime-JIT candidate restores `build-essential`, uses
explicit GCC/G++ paths, and measures 14.97 GB. It passes package,
upstream-import, real-checkpoint, wrapper, static workflow, and native Triton
extension-build checks, but has not been published because cold-cache Docker
and Singularity GPU/JIT tests are still required. Both modules now select
`phenogpt2_0.4`.

```bash
# The persistent runtime cache is optional, but must already exist when supplied.
mkdir -p /data/phenogpt2-cache/image-0.3_models-v1_h100
nextflow run main.nf -profile local_docker \
  --bam sample.bam --ref_fa /refs/hg38.fa --note note.txt \
  --phenotype_extractor phenogpt2 --GPU yes \
  --phenogpt2_negation yes \
  --phenogpt2_model_host_path /data/models/phenogpt2-v1/new_model \
  --phenogpt2_negation_model_host_path /data/models/phenogpt2-negation \
  --phenogpt2_embedding_model_host_path /data/models/phenogpt2-embedding \
  --phenogpt2_cache_host_path /data/phenogpt2-cache/image-0.3_models-v1_h100
```

Each path may be either a complete standalone model directory or a complete
Hugging Face cache repository containing `refs/main` and its referenced
snapshot. No checksum-manifest filename or model filename is imposed by
PipeVar. PipeVar preflight checks only that the configured host directories
exist. Inside the container, the upstream Transformers and SentenceTransformer
loaders determine whether each mounted model is usable. For offline embedding
loading, the wrapper discovers the repository identifier used by the pinned
upstream source and creates the corresponding runtime-cache alias; it does not
duplicate a model identifier. Upstream prompts, thresholds, GPU placement, and
embedding-based HPO verification remain unchanged.

The runtime cache also contains CUDA, TorchInductor, Triton, and vLLM
subdirectories. Keep it writable; use a persistent cache only when it is
namespaced for the image, checkpoint version, and GPU class.

### Setup Script

ANNOVAR must be downloaded through the ANNOVAR registration process:

- https://www.openbioinformatics.org/annovar/annovar_download_form.php

Then run:

```bash
# Full setup
./setup.sh

# Light PhenoSV setup
./setup.sh light
```

By default, setup expects ANNOVAR at `./annovar` and PhenoSV resources under `./PhenoSV_model`. You can override both locations:

```bash
./setup.sh --annovar-dir=/shared/apps/annovar --phenosv-dir=/shared/data/PhenoSV_model
```

The setup script prepares required assets and can persist the default execution profile plus bind source paths in `nextflow.config`.

Non-interactive setup example:

```bash
./setup.sh --non-interactive --profile=local_docker \
  --annovar-dir=/data/annovar \
  --phenosv-dir=/data/PhenoSV_model \
  --annovar-bind=/data/annovar \
  --phenosv-bind=/data/PhenoSV_model
```

## Input Modes

### Single BAM/CRAM

Required:

- `--bam <FILE>`
- `--ref_fa <FILE>`
- one phenotype source:
  - `--note <FILE>` for a clinical note processed by PhenoTagger
  - `--hpo <FILE>` for an HPO term file

Optional:

- `--mode snp` or `--mode sv` to run only one branch.
- omit `--mode` to run the full supported branch set for the selected sequencing type.
- `--type ont`, `--type pacbio`, or `--type short` to select the sequencing path.

BAM inputs require a `.bai` index. CRAM inputs require a `.crai` index and the matching reference FASTA.

### Single VCF

Required:

- `--vcf <FILE>`
- `--ref_fa <FILE>`
- `--mode snp` or `--mode sv`
- one phenotype source: `--note <FILE>` or `--hpo <FILE>`

Single VCF mode re-annotates and prioritizes an existing VCF. Mitochondrial analysis and xTEA require BAM/CRAM input and are not available in VCF-only mode.

### Annotated SNV And Annotated SV

Annotated SNV mode starts from ANNOVAR multianno SNV outputs.

Single-sample annotated SNV required:

- `--annotated_snv yes`
- `--annovar_txt <ANNOVAR multianno TXT>`
- `--vcf <matching ANNOVAR multianno VCF>`
- one phenotype source: `--note <FILE>` or `--hpo <FILE>`

To reuse annotated SNVs while calling SV/STR and optional mito from short-read BAM/CRAM:

```bash
nextflow run main.nf \
  --annotated_snv yes \
  --annovar_txt sample.hg38_multianno.txt \
  --vcf sample.hg38_multianno.vcf \
  --bam sample.bam \
  --ref_fa ref.fa \
  --type short \
  --hpo sample.hpo.txt
```

To import annotated SVs instead of calling SVs from BAM/CRAM, add:

- `--annotated_sv yes`
- `--annovar_sv_vcf <SV multianno VCF>`

An alignment is optional when both SNV and SV inputs are already annotated:

```bash
nextflow run main.nf \
  --annotated_snv yes \
  --annotated_sv yes \
  --annovar_txt sample.hg38_multianno.txt \
  --vcf sample.hg38_multianno.vcf \
  --annovar_sv_vcf sample.sv.hg38_multianno.vcf \
  --hpo sample.hpo.txt
```

This VCF-only route performs combined SNV/SV prioritization and writes the final HTML report, but it does not run STR or mitochondrial analysis and does not require `--ref_fa`.

Annotated SNV mode is SNP-led. Use `--mode snp` or omit `--mode`; `--mode sv` is not supported for this input family.

### Legacy CSV

Legacy CSV mode uses `sample,file_path,note_path` rows.

BAM/CRAM CSV required:

- `--input_csv <FILE>`
- `--bam true`
- `--ref_fa <FILE>`

VCF CSV required:

- `--input_csv <FILE>`
- `--vcf true`
- `--ref_fa <FILE>`
- `--mode snp` or `--mode sv`

Shared legacy CSV columns:

| Column | Required | Description |
| --- | --- | --- |
| `sample` | yes | Output/sample identifier |
| `file_path` | yes | BAM/CRAM or VCF path |
| `note_path` | proband only with de novo; otherwise yes | Clinical note by default, or HPO file when `--note no` |
| `age_of_onset` | no | Per-sample age for prioritization |
| `age` | no | Alternate age column used only when `age_of_onset` is absent |
| `sex` | no | Per-sample sex metadata; override name with `--sex_column` |

Age handling:

- `age_of_onset` is preferred if both age columns are present.
- Empty age is allowed and treated as not provided.
- Non-empty age must be `<integer>`, `<integer>d`, `<integer>m`, or `<integer>y`.
- Integer-only ages are normalized to years, for example `7` becomes `7y`.

Sex handling:

- Default sex column name is `sex`; override with `--sex_column <STRING>`.
- Values are normalized to lowercase and must be `unknown`, `male`, or `female`.
- Empty or missing values are treated as `unknown`.

Phenotype handling:

- By default, `note_path` is treated as a clinical note and PhenoTagger runs.
- With `--phenotype_extractor phenogpt2 --GPU yes`, clinical notes are processed with GPU-backed PhenoGPT2 instead of PhenoTagger.
- PhenoGPT2 requires `--phenogpt2_model_host_path /absolute/versioned/new_model`; with `--phenogpt2_negation yes`, it also requires the negation and embedding model host paths.
- With `--note no`, `note_path` is treated as an HPO file and PhenoTagger is skipped.

### Unified CSV

Unified CSV mode uses `input_kind` and does not combine `--input_csv` with legacy `--bam true` or `--vcf true` flags.

Required base columns:

| Column | Required | Values |
| --- | --- | --- |
| `sample` | yes | Unique sample identifier |
| `input_kind` | yes | `annotated_snv`, `vcf_snv`, `vcf_sv`, `bam_ngs`, `cram_ngs` |
| `phenotype_path` | yes | Clinical note or HPO file path |
| `phenotype_format` | yes | `clinical_note` or `hpo` |

Optional/conditional columns:

| Column | Used by | Description |
| --- | --- | --- |
| `age_of_onset` | all CSV prioritization flows | Per-sample onset age |
| `sex` | all CSV prioritization flows | Per-sample sex metadata |
| `snv_txt_path` | `annotated_snv` | ANNOVAR multianno TXT |
| `snv_vcf_path` | `annotated_snv` | Matching ANNOVAR multianno VCF |
| `sv_vcf_path` | `annotated_snv` combined SNV/SV modes | Annotated SV VCF; may be used with or without alignment input |
| `vcf_path` | `vcf_snv`, `vcf_sv` | Existing VCF input |
| `alignment_path` | `bam_ngs`, `cram_ngs`, annotated hybrid modes | BAM/CRAM input |
| `alignment_index_path` | alignment-backed modes | Optional explicit BAM/CRAM index path |

Current unified manifest constraints:

- Samples must be unique.
- Mixed `input_kind` sets are limited to supported homogeneous groups, except mixed `bam_ngs` and `cram_ngs`.
- Annotated SNV rows either all provide `alignment_path` or none do.
- Annotated SNV rows either all provide `sv_vcf_path` or all leave it blank.
- Mixed `phenotype_format` values are not supported for non-annotated CSV modes.

The four homogeneous `annotated_snv` manifest combinations are SNV only, annotated SNV+SV without alignment, annotated SNV plus alignment with SV calling, and annotated SNV+SV plus alignment. The no-alignment combinations cannot run STR or mitochondrial analysis.

Run `scripts/generate_input_csv.sh` to build legacy or unified CSV manifests interactively. Its annotated-SNV update action can add `sv_vcf_path` to manifests with or without alignments; choose a separate output path so the source CSV remains unchanged.

## Feature Options

### Mitochondrial Analysis

Enable with `--mito yes`.

- Supported only for BAM/CRAM input.
- Supported with `--mode snp` or with `--mode` omitted.
- Not supported with `--mode sv`.
- Short reads use Mutect2.
- Long reads use a mito-specific Clair3 call step followed by VCF postprocessing.
- Long-read mito is unavailable with `--light yes` because that path selects NanoCaller.

The mito branch emits separate mito outputs and does not modify nuclear `.prio.vcf` outputs. Bundled mito evidence sources include MITOMAP, MitoTip, t-APOGEE, and MitImpact.

### xTEA Mobile-Element Analysis

Enable with `--xtea yes`.

- Supported only for `--type short`.
- Supported only for BAM/CRAM input.
- Supported with `--mode sv` or with `--mode` omitted.
- Not supported with `--mode snp`.
- Runs per sample in single-sample and CSV batch BAM/CRAM modes.

The xTEA image is expected to contain the `xtea` command, xTEA scripts under `/opt/xTea/xtea`, the repeat library under `/opt/xtea/rep_lib_annotation`, and the GENCODE GFF3 annotation at `/opt/xtea/gencode.gff3`.

### CNVnator

CNVnator is enabled by default for short-read SV/all-NGS calling with `--cnvnator yes`.

Long-read SV discovery uses Sniffles. Full ONT and PacBio workflows pass the complete
Sniffles VCF, including supporting read names, to LongPhase. Annotation, common-SV,
de novo, and phenotype filtering remain a separate evidence branch. PipeVar uses the
whole-event PhenoSV score from the `Elements=SV` row, obtains gene symbols from the
curated ANNOVAR `Gene.refGene` annotation, and joins both to LongPhase's phased
Sniffles output by exact VCF ID.

At the SURVIVOR/PhenoSV boundary, PipeVar_mito creates four staged files: a canonical
simple-SV BED with `DEL`, `DUP`, `INV`, and `INS`; a temporary PhenoSV BED using
`deletion`, `duplication`, `inversion`, and `insertion`; a headerless nine-column BEDPE
for BND/TRA adjacencies; and a member TSV mapping each scored adjacency to its original
VCF ID or reciprocal pair of IDs. Simple BED and BEDPE partitions are scored separately
and merged beneath one output header. Reciprocal BND scores are expanded back to both
original IDs, so downstream VCF reconstruction retains the original ALT, `SVTYPE=BND`,
`MATEID`, genotype, and record IDs. Valid BNDs without `MATEID` are scored as singleton
adjacencies.

PhenoSV-light continues to score BEDPE input but emits a reduced-translocation-accuracy
warning. The pinned upstream PhenoSV checkout is not modified: PipeVar sends
`duplication` at the tool boundary and restores canonical `DUP` in its own output, while
the pinned upstream implementation still applies its deletion-like internal duplication
transformation. The reserved `survivor_0.1` and `phenosv_0.2` images must be built before
these paths can run.

This repair also reserves `rankscore_0.2.22`, `longphase_0.2.33`,
`validate_preannotated_annovar_pair_0.3`, and `mito_annotation_0.4.2`. These tags are
intentionally unpublished in this source change;
a post-build mixed DEL/DUP/INV/INS/BND smoke test is a release gate before publication.

### ClinVar And Compound-Heterozygous Semantics

ClinVar filtering reads the `CLNSIG` column by header name rather than searching whole
ANNOVAR rows. Accepted values are canonicalized to `Pathogenic`, `Likely_pathogenic`, or
`Pathogenic/Likely_pathogenic`; low-penetrance P/LP assertions collapse to the matching
base class. Conflicting, VUS, benign, and risk-only assertions are excluded. The
`--include_clinvar_report` toggle controls only ClinVar-exclusive report entries and does
not alter P/LP classification or upgrade likely pathogenic calls to pathogenic.

Compound-heterozygous evaluation follows this phase contract:

| Pair state | Result |
| --- | --- |
| Same valid PS, opposite phased orientations (`1|0` + `0|1`) | Confirmed trans; accepted |
| Same valid PS, same phased orientation | Proven cis; always rejected |
| Different valid PS values | Unresolved; requires `--allow_unphased_comphet yes` |
| Either PS missing or malformed, opposite phased orientations | GT-orientation trans; accepted |
| Either PS missing or malformed, same phased orientation | GT-orientation cis; always rejected |
| Slash-unphased pair | Unresolved; requires the opt-in |
| Mixed phased and slash-unphased pair | Unresolved; requires the opt-in |

Assignment resolves `GT` and optional `PS` through FORMAT field names, including
`GT:PS` and `PS:GT`. Valid integer PS values and original VCF IDs are preserved;
missing PS is accepted, while malformed PS is warned, omitted, and evaluated using GT.
GT orientation without a shared phase block is a heuristic: orientation may flip between
blocks, so GT-fallback trans/cis labels are not equivalent to shared-PS confirmation.

### Common-SV Filtering

Common-SV filtering removes common SVs before PhenoSV and final SV prioritization. It is enabled by default in `nextflow.config` with `--common_sv_filter yes`.

Matching thresholds are controlled by `--common_sv_af`, `--common_sv_reciprocal_overlap`, `--common_sv_distance`, `--common_sv_ins_distance`, and `--common_sv_ins_identity`.

### De Novo Filtering

CSV de novo filtering is disabled by default with `--denovo_filter no`.

When enabled, CSV family metadata is used to compare cohort VCF calls before any sample-specific Phen2Gene or target-region restriction. Called SNV and SV records are filtered before ANNOVAR; imported ANNOVAR SNV TXT/VCF pairs are filtered together. Every CSV row must have a `proband`, `father`, `mother`, or `sibling` role; each family must contain one proband and at least one parent. Sample identifiers must be globally unique and use only letters, numbers, `.`, `_`, and `-`.

Parents are variant controls only: their phenotype, age, and sex cells may be blank, and only probands continue through phenotype extraction, sex-aware prioritization, repeat/mitochondrial analysis, and final reporting. With both de novo and sex metadata enabled, inherited maternal chrX and paternal chrY calls are removed first; retained male chrX/XLR and chrY/AR calls can then receive hemizygous prioritization. Female and `unknown` probands retain the existing non-hemizygous behavior.

Example combined columns:

```csv
sample,file_path,note_path,family_id,role,vcf_sample,sex
child,child.bam,child.hpo.txt,F1,proband,CHILD,male
mother,mother.bam,,F1,mother,MOTHER,
father,father.bam,,F1,father,FATHER,
```

Configure role, family, and sample mapping columns with `--denovo_role_column`, `--denovo_family_column`, and `--denovo_vcf_sample_column`. Published binding manifests preserve sample-to-file identity, while `denovo.snv.summary.tsv` and `denovo.sv.summary.tsv` report the parents used, parental variants loaded, proband variants examined, and kept/removed counts. A zero-match family produces a warning rather than failing because it can be biologically valid.

### Light Mode

`--light yes` switches caller/model choices internally:

| Branch | Default | `--light yes` |
| --- | --- | --- |
| short-read SNP | DeepVariant | GATK HaplotypeCaller |
| long-read SNP | Clair3 | NanoCaller |
| SV prioritization | PhenoSV default model | PhenoSV-light config |

### Targeted And Gene Filtering

`--target yes` restricts SNP calling to phenotype-derived gene regions where supported.

`--phen2gene_filter <INT>` sets the number of top Phen2Gene genes retained for targeted mode.

SV annotation is not restricted by `--target` or Phen2Gene scores. Called and imported
SVs are retained at this stage when ANNOVAR reports `Func.refGene=exonic`; optional
common-SV filtering and PhenoSV scoring still apply downstream.

`--gene <SYMBOLS|FILE>` restricts final prioritization to comma-separated gene symbols or a one-gene-per-line file.

## Parameters

Defaults below come from `nextflow.config`.

### Input, Output, And Runtime

| Parameter | Default | Description |
| --- | --- | --- |
| `--bam` | `null` | Single BAM/CRAM input |
| `--vcf` | `null` | Single VCF input |
| `--input_csv` | `null` | CSV manifest |
| `--annovar_txt` | `null` | Single-sample ANNOVAR multianno TXT for annotated SNV mode |
| `--annovar_sv_vcf` | `null` | Single-sample ANNOVAR SV multianno VCF for annotated SV import |
| `--ref_fa` | `null` | Reference FASTA |
| `--out_prefix` | `PipeVar` | Single-sample output prefix |
| `--output_directory` | launch directory | Publish directory |
| `--mode` | `null` | `snp` or `sv`; omit for all supported branches |
| `--type` | `ont` | `ont`, `pacbio`, or `short` |
| `--light` | `null` | Use lightweight callers/models when set to `yes` |
| `--genome` | `hg38` | Genome build for default ExpansionHunter catalog |
| `--expansionhunter_variant_catalog` | `null` | Optional ExpansionHunter catalog override |
| `--sex_column` | `sex` | Optional CSV sex metadata column name |
| `--annovar_host_path` | `${projectDir}/annovar` | Host ANNOVAR path mounted into containers |
| `--phenosv_host_path` | `${projectDir}/PhenoSV_model` | Host PhenoSV resource path mounted into containers |
| `--help` | `false` | Print compact help and exit |

### Phenotype And Prioritization

| Parameter | Default | Description |
| --- | --- | --- |
| `--note` | `null` | Clinical note path, or `no` in CSV mode to treat `note_path` as HPO |
| `--hpo` | `null` | HPO term file |
| `--gnomad` | `0.0001` | Maximum gnomAD AF for SNP prioritization |
| `--inheritance_mode` | `ml` | `ml`, `omim`, or `gnomad`; `gnomad` maps to LOEUF fallback lists |
| `--include_clinvar_report` | `yes` | Include ClinVar-only calls in final prioritized outputs |
| `--allow_unphased_comphet` | `no` | Allow unresolved AR pairs with different valid PS, slash-unphased GT, or mixed phased/unphased GT; GT fallback handles missing/malformed PS |
| `--prioritize_sv_only` | `no` | In combined final prioritization, report only SV/PhenoSV evidence |
| `--rankscore` | `0.50` | Minimum RankScore cutoff |
| `--rankscore_softwares` | `null` | Comma-separated RankScore software list; null means all built-in tools |
| `--rankvar` | `0.05` | Minimum RankVar score |
| `--phenosv_score` | `0.50` | Minimum PhenoSV score |
| `--gq` | `20` | Minimum genotype quality |
| `--ad` | `15` | Minimum allele depth |
| `--nanocaller_dp` | `20` | Minimum NanoCaller depth used by RankVar only when GQ is missing; this is a coverage threshold, not a GQ equivalent |
| `--phen2gene_filter` | `500` | Number of top Phen2Gene genes for targeted mode |
| `--gene` | `null` | Comma-separated genes or one-gene-per-line file |
| `--target` | `null` | Enable phenotype-derived targeted calling with `yes` |
| `--phenotype_extractor` | `phenotagger` | Clinical-note extractor: `phenotagger` or GPU-backed `phenogpt2` |
| `--GPU` | `no` | Shared GPU mode for DeepVariant GPU and PhenoGPT2 |
| `--gpu_backend` | `singularity` | GPU container backend: `singularity` uses `--nv`, `docker` uses `--gpus 1` |
| `--gpu_cpus` | `16` | CPU threads assigned to GPU-backed DeepVariant/PhenoGPT2 processes |
| `--gpu_cluster_options` | `--gres=gpu:1` | Scheduler options applied to GPU-backed processes |
| `--phenogpt2_batch_size` | `1` | PhenoGPT2 inference batch size; currently required to remain `1` |
| `--phenogpt2_chunk_batch_size` | `1` | PhenoGPT2 chunk batch size; currently required to remain `1` |
| `--phenogpt2_wc` | `0` | PhenoGPT2 word-count chunking; `0` disables chunking |
| `--phenogpt2_attn_implementation` | `eager` | PhenoGPT2 attention implementation |
| `--phenogpt2_negation` | `no` | Enable the full upstream negation and embedding-based HPO verification pass |
| `--phenogpt2_model_host_path` | `null` | Required path to an existing `new_model` directory; mounted read-only and validated by the container wrapper |
| `--phenogpt2_negation_model_host_path` | `null` | Required with negation; existing negation-model directory mounted read-only and validated by the wrapper |
| `--phenogpt2_embedding_model_host_path` | `null` | Required with negation; existing embedding-model directory mounted read-only and validated by the wrapper |
| `--phenogpt2_cache_host_path` | `null` | Optional absolute canonical path to a pre-created writable persistent cache; otherwise each task uses its own work-directory cache |
| `--phenogpt2_max_forks` | `1` | Maximum concurrent PhenoGPT2 tasks; currently required to remain `1` |

For NanoCaller-selected long-read paths, RankVar uses genuine GQ when present and
does not rescue a low GQ with DP. If GQ is missing, DP must be at least
`--nanocaller_dp` and first-ALT AD must be greater than `--ad`. A passing GQ can
retain NanoCaller indels without AD. Other callers do not receive the DP fallback.

### Caller And Feature Toggles

| Parameter | Default | Description |
| --- | --- | --- |
| `--annotated_snv` | `no` | Use pre-annotated ANNOVAR TXT + VCF SNV input |
| `--annotated_sv` | `no` | Use pre-annotated ANNOVAR SV VCF with annotated SNV + short-read BAM/CRAM |
| `--cnvnator` | `yes` | Add CNVnator to short-read SV/all-NGS calling |
| `--cnvnator_bin_size` | `100` | CNVnator read-depth bin size |
| `--xtea` | `no` | Add xTEA mobile-element calling to short-read SV/all-NGS BAM/CRAM paths |
| `--mito` | `no` | Add mitochondrial analysis for BAM/CRAM input |

### Filtering And Reporting

| Parameter | Default | Description |
| --- | --- | --- |
| `--common_sv_filter` | `yes` | Remove common SVs before PhenoSV and final SV prioritization |
| `--common_sv_af` | `0.01` | Minimum AF used to treat baked common SV database records as common |
| `--common_sv_reciprocal_overlap` | `0.5` | Minimum reciprocal overlap for interval SV matching |
| `--common_sv_distance` | `1000` | Breakpoint fallback distance for interval SV matching |
| `--common_sv_ins_distance` | `500` | Insertion position window for common matching |
| `--common_sv_ins_identity` | `0.5` | Insertion sequence identity threshold when inserted sequence is available |
| `--denovo_filter` | `no` | Filter proband SNV/SV calls against parents before phenotype-specific restriction |
| `--denovo_role_column` | `role` | CSV role column for proband/father/mother/sibling |
| `--denovo_family_column` | `family_id` | CSV family grouping column |
| `--denovo_vcf_sample_column` | `vcf_sample` | Optional CSV column for VCF sample names |
| `--denovo_sv_min_reciprocal_overlap` | `0.50` | SV parent/proband reciprocal-overlap threshold |
| `--denovo_exclude_contigs` | `MT,M,chrM,chrMT` | Contigs excluded from de novo filtering |
| `--mito_contig` | `chrM` | Preferred mitochondrial contig alias |
| `--mito_min_vaf` | `0.01` | Mito prioritization VAF floor |
| `--mito_min_depth` | `50` | Mito prioritization depth floor |
| `--mito_min_alt_reads` | `5` | Mito prioritization alternate-read floor |
| `--mito_gui_min_af` | `0.5` | Strict mtDNA allele fraction for final report mitochondrial rows |
| `--mito_gui_min_apogee2` | `0.5` | Strict APOGEE2 score for final report mitochondrial rows |
| `--mito_gui_min_mitotip` | `12.66` | Strict MitoTip score for final report mitochondrial rows |

## Example Commands

### Single-Sample Long-Read Full Analysis

```bash
nextflow run main.nf \
  -profile standard \
  --bam /data/p1.bam \
  --ref_fa /refs/hg38.fa \
  --note /data/p1_note.txt \
  --out_prefix p1 \
  --type ont
```

### Single-Sample Short-Read Full Analysis With Mito

```bash
nextflow run main.nf \
  -profile standard \
  --bam /data/p2.bam \
  --ref_fa /refs/hg38.fa \
  --hpo /data/p2_hpo.txt \
  --out_prefix p2 \
  --type short \
  --mito yes
```

### Single-Sample Short-Read SV Analysis With xTEA

```bash
nextflow run main.nf \
  -profile standard \
  --bam /data/p2.bam \
  --ref_fa /refs/hg38.fa \
  --hpo /data/p2_hpo.txt \
  --out_prefix p2_sv \
  --type short \
  --mode sv \
  --xtea yes
```

### Single-Sample VCF SNP Re-Annotation And Prioritization

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

### CSV Batch Short-Read SV Mode With xTEA

```bash
nextflow run main.nf \
  -profile slurm_singularity \
  --input_csv /data/samples.csv \
  --bam true \
  --note no \
  --ref_fa /refs/hg38.fa \
  --type short \
  --mode sv \
  --xtea yes
```

### CSV Batch VCF Mode

```bash
nextflow run main.nf \
  -profile local_singularity \
  --input_csv /data/sv_samples.csv \
  --vcf true \
  --mode sv \
  --ref_fa /refs/hg38.fa
```

## Outputs

Outputs are published to `--output_directory`. Exact files depend on `--mode`, `--type`, input kind, and enabled feature toggles.

### SNP-Related Outputs

- Caller outputs:
  - `*.deepvariant.vcf.gz` for short-read default SNP calling.
  - `*.recal.vcf.gz` for the short-read HaplotypeCaller/light path.
  - `*.clair3.vcf.gz` for long-read default SNP calling.
  - `*.nanocaller.vcf.gz` for the long-read light path.
  - `*.mito.vcf.gz` for the mito branch.
- Annotation and prioritization outputs:
  - `*.clinvar.txt`
  - `*.rank_var.tsv`
  - `*.rankscore_filtered.tsv`
  - ANNOVAR intermediate/final files such as `*.hg38_multianno.*`
  - `*.mito.annotated.tsv`
  - `*.mito.annotated.vcf.gz`
  - `*.mito.prioritized.tsv`

### SV-Related Outputs

- Short-read SV/MEI:
  - `*_manta.vcf`
  - `*_xtea.vcf` when `--xtea yes`
  - `*.shortread_sv.merged.vcf`
  - `*.shortread_sv.truvari_collapsed.vcf`
- Long-read SV:
  - `*.sniffles.vcf` (includes `RNAMES` for LongPhase phasing)
- Downstream SV prioritization:
  - `*.exonic.vcf`
  - `*.phenosv.filtered.tsv`, with one whole-SV row per passing event and columns
    `SV_ID`, `CHROM`, `START`, `END`, `SVTYPE`, `PATHOGENICITY`, `PHEN2GENE`,
    `PHENOSV_SCORE`, and `PHENOSV_TYPE`
  - Final prioritized VCFs use ANNOVAR genes with the event-level `PHENO_SCORE`;
    PhenoSV gene-row scores are not used as independent evidence.

### Repeat Expansion Outputs

- Short-read ExpansionHunter outputs:
  - `*.json`
  - `*.eh.tsv`
- Long-read NanoRepeat outputs:
  - `*_nanoRepeat_output.tsv`
  - related summary files

### Phenotype Outputs

- `*_phenotagger_patient_hpo.txt`
- `*_phenogpt2_patient_hpo.txt` when `--phenotype_extractor phenogpt2`
- Phen2Gene ranking outputs such as `*_phen2gene*`

## Notes And Pitfalls

- `--input_csv` legacy mode requires either `--bam true` or `--vcf true`.
- Unified CSV mode should not be combined with `--bam true` or `--vcf true`.
- Single-file mode requires `--note <FILE>` or `--hpo <FILE>`.
- Single VCF mode requires `--mode snp` or `--mode sv`.
- Reference index `${ref}.fai` must exist when `--ref_fa` is supplied.
- BAM/CRAM indexes must exist for alignment-driven paths.
- Short-read mitochondrial analysis also requires a reference dictionary and BWA sidecars: `.dict`, `.amb`, `.ann`, `.bwt`, `.pac`, and `.sa`.
- The DRAGEN CRAM compatibility path uses `RevertSam --RESTORE_HARDCLIPS false`.
- xTEA is intended for short-read WGS MEI discovery/genotyping and requires indexed BAM/CRAM input.
- Long-read raw Sniffles VCFs contain supporting read identifiers and can be larger than VCFs produced without `--output-rnames`.
- Compound-heterozygous calls use matching valid PS when available and otherwise fall back to pipe-phased GT orientation. Different valid PS values, slash-unphased GT, and mixed phased/unphased GT require `--allow_unphased_comphet yes`; same-orientation cis pairs are rejected unless their valid PS values differ, in which case the pair is unresolved.
- If using Singularity/Docker profiles, ensure `--annovar_host_path` and `--phenosv_host_path` point to valid host locations.
- External-model PhenoGPT2 support is text-only with optional full upstream negation and `--phenogpt2_wc 0`; vision, training, and BERT chunk filtering are not provisioned.
- Keep all mounted PhenoGPT2 checkpoints immutable and use a new versioned directory for every model change. PipeVar checks that each configured host directory exists, and the container wrapper uses the upstream loaders to validate model usability. Because Nextflow does not hash model contents during preflight, use a fresh work directory or `-resume` only when the mounted model versions are unchanged.
- On HPC, model and persistent-cache paths must exist at the identical absolute path on the submit node and every GPU compute node. Persistent caches must be trusted, writable, and namespaced by PhenoGPT2 image, model version, and GPU class.
- Negation validation targets one scheduler-isolated 80 GB GPU. A 40 GB GPU is supported only after cold-cache and long-note testing completes without OOM. External mounting reduces image transfer size, not checkpoint distribution size or GPU-memory demand.

## Additional Documentation

- `docs/USAGE.md` for additional usage notes.
- `docs/WORKFLOW.md` for workflow-level documentation.
- `docs/PIPEVAR_NUCLEAR_FIGURES.md` for focused nuclear workflow figures.
- `scripts/generate_input_csv.sh` for generating legacy or unified CSV manifests.

## Software Components

### SNP Calling

- DeepVariant
- GATK HaplotypeCaller
- Clair3
- NanoCaller

### SV/CNV Calling And Prioritization

- Sniffles
- Manta
- CNVnator
- xTEA
- SURVIVOR
- Truvari
- PhenoSV
- ANNOVAR SV annotation

### Repeat Expansion

- ExpansionHunter
- NanoRepeat

### Annotation, Ranking, And Phenotype

- ANNOVAR
- RankVar
- RankScore filtering
- Phen2Gene
- PhenoTagger
- Longphase prioritization helpers

## Status

PipeVar_mito is under active development. `nextflow.config` defines runtime defaults, `README.md` and linked docs describe supported user-facing modes, and `nextflow run main.nf --help` provides a compact command-line quick reference.
