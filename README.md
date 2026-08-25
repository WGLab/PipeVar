# PipeVar

PipeVar analyzes short- or long-read sequencing data and ranks DNA changes that
may help explain a patient's symptoms. It accepts aligned reads, existing
variant calls, or prepared annotations and produces prioritized candidates for
review.

PipeVar 0.5.0 supports:

- phenotype-guided rare-disease analysis;
- BAM/CRAM, VCF, and prepared ANNOVAR inputs;
- short- and long-read nuclear analysis, repeat analysis, and optional
  mitochondrial analysis; and
- single-sample or batch execution with prioritized files and reports.

## Documentation

| I want to... | Read... |
| --- | --- |
| Install PipeVar or configure containers | [Installation and setup](docs/INSTALLATION.md) |
| Prepare a reference, sample, phenotype file, or sample sheet | [Input guide](docs/INPUTS.md) |
| Choose a workflow or find more run examples | [Running PipeVar](docs/USAGE.md) |
| Look up an option and its default | [Parameter reference](docs/PARAMETERS.md) |
| Understand the analysis stages | [Workflow guide](docs/WORKFLOW.md) |
| Find and interpret published files | [Output guide](docs/OUTPUTS.md) |

## Choose a workflow

In this documentation, BAM and CRAM are aligned-read files, VCF contains
existing variant calls, and ANNOVAR multianno files are annotation tables
created by ANNOVAR. **Small variants** are DNA-letter changes and short
insertions/deletions; **structural/copy-number variants** are larger
rearrangements, gains, or losses. **Combined** means analyzing both groups.

| Starting point | Sequencing | Available analysis | Execution |
| --- | --- | --- | --- |
| BAM/CRAM | Short read, Oxford Nanopore, or PacBio | Small variants, structural/copy-number variants, or combined | Single sample or sample sheet |
| VCF | Existing small-variant or structural-variant calls | Re-annotation and prioritization | Single sample or sample sheet |
| ANNOVAR multianno files | Prepared small-variant annotations | Re-prioritization, with supported batch combinations | Single sample or typed sample sheet |

See [Running PipeVar](docs/USAGE.md) for route restrictions and
[Input guide](docs/INPUTS.md) for required files and indexes.

## Quick start

### 1. Provide the required software and data

Install Nextflow with Java and either Docker or Singularity. Obtain ANNOVAR
through its registration process, and prepare an hg38 reference for your data.
Aligned-read examples also require a discoverable BAM/CRAM index and
`<reference>.fai`; see [reference files](docs/INPUTS.md#reference-files) and
[alignment indexes](docs/INPUTS.md#alignment-indexes).

Clone PipeVar and enter the repository:

```bash
git clone https://github.com/WGLab/PipeVar.git PipeVar
cd PipeVar
```

### 2. Prepare PipeVar resources

`setup.sh` is PipeVar's setup assistant. It validates an existing licensed
ANNOVAR installation, downloads the hg38 ANNOVAR databases used by PipeVar,
downloads PhenoSV resources and the Phen2Gene knowledge base, adjusts the
resource paths, and configures `profiles.standard` for the selected container
backend. Run it from the repository root.

The script does not install Nextflow, Java, Docker, Singularity, or ANNOVAR
itself. It also does not provide a reference FASTA and indexes, PhenoGPT2
models, a custom repeat catalog, or optional HmtVar data.

For guided setup:

```bash
bash setup.sh
```

For reproducible local Docker setup:

```bash
bash setup.sh --non-interactive \
  --profile=local_docker \
  --annovar-dir=/absolute/path/to/annovar \
  --phenosv-dir=/absolute/path/to/PhenoSV_model
```

For reproducible local Singularity setup:

```bash
bash setup.sh --non-interactive \
  --profile=local_singularity \
  --annovar-dir=/absolute/path/to/annovar \
  --phenosv-dir=/absolute/path/to/PhenoSV_model
```

Setup rewrites resource paths and `profiles.standard`. Review the detailed
[installation guide](docs/INSTALLATION.md) before rerunning it or configuring
separate bind paths, GPU execution, or a cluster.

The setup `--profile=...` option chooses the backend written into
`profiles.standard`, which is used when a run does not specify a Nextflow
profile. The examples below deliberately use `-profile local_docker` or
`-profile local_singularity` to select an existing named profile directly.

### 3. Run one sample with Docker

Replace `/data`, `/refs`, and `/results` with local absolute paths, and ensure
the output location is writable.

```bash
nextflow run main.nf \
  -profile local_docker \
  --bam /data/proband.bam \
  --ref_fa /refs/hg38.fa \
  --hpo /data/proband.hpo.txt \
  --type short \
  --out_prefix proband \
  --output_directory /results/proband
```

| Argument | Meaning |
| --- | --- |
| `-profile local_docker` | Nextflow execution profile; run tasks locally in Docker containers |
| `--bam` | PipeVar input parameter for the BAM or CRAM aligned-read file |
| `--ref_fa` | Matching hg38 reference FASTA |
| `--hpo` | File containing standardized Human Phenotype Ontology codes for the patient's clinical findings |
| `--type short` | Select the short-read workflow |
| Omitted `--mode` | Run the supported combined nuclear analysis |
| `--out_prefix` | Prefix for this sample's output files |
| `--output_directory` | Directory that receives published results |

Nextflow options use one hyphen, such as `-profile`. PipeVar parameters use two
hyphens, such as `--bam` and `--ref_fa`.

### 4. Run the same sample with Singularity

Run the same command with `-profile local_singularity` in place of
`-profile local_docker`. All PipeVar input parameters keep the same meaning.

## Results and next steps

The example publishes ranked candidate variants and genes and, for supported
combined routes, an HTML report. Open the HTML report first for a human-readable
summary. These are candidates for expert review rather than a diagnosis. See
the [output guide](docs/OUTPUTS.md) for exact filenames and conditional results.

For batch runs, VCF input, long-read analysis, prepared annotations,
mitochondrial analysis, and cluster examples, continue with
[Running PipeVar](docs/USAGE.md).

## Status and support

PipeVar 0.5.0 is under active development. This README and the linked guides
define the documented operator interface.

Report reproducible problems through the
[PipeVar GitHub issue tracker](https://github.com/WGLab/PipeVar/issues). Include
the PipeVar revision, command, profile, Nextflow version, container backend,
reference build, and relevant report or trace excerpts.
