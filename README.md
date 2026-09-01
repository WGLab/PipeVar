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
| Plan compute capacity or review the 30× example | [Resources and benchmark](docs/BENCHMARKS.md) |
| Prepare a reference, sample, phenotype file, or sample sheet | [Input guide](docs/INPUTS.md) |
| Choose a workflow or find more run examples | [Running PipeVar](docs/USAGE.md) |
| Look up an option and its default | [Parameter reference](docs/PARAMETERS.md) |
| Understand the analysis stages | [Workflow guide](docs/WORKFLOW.md) and [nuclear figures](docs/PIPEVAR_NUCLEAR_FIGURES.md) |
| Find and interpret published files | [Output guide](docs/OUTPUTS.md) |

## Choose a workflow

In this documentation, BAM and CRAM are aligned-read files, VCF contains
existing variant calls, and ANNOVAR multianno files are annotation tables
created by ANNOVAR. **Small variants** are single-letter changes and short
insertions/deletions; **structural/copy-number variants** are larger
rearrangements, gains, or losses.

For aligned reads, make two explicit choices:

1. Set the sequencing technology. Do not rely on the default.

   | Your data | Type this |
   | --- | --- |
   | Illumina or other short reads, including whole-exome (WES) or whole-genome sequencing (WGS) | `--type short` |
   | Oxford Nanopore | `--type ont` |
   | PacBio | `--type pacbio` |

2. Select the analysis.

   | What you want | Type this |
   | --- | --- |
   | Small variants only | `--mode snp` |
   | Structural/copy-number variants and repeats | `--mode sv` |
   | Combined small and structural/copy-number analysis | Do not include `--mode` |

3. Select one sample or a batch.

   | Execution | Type this |
   | --- | --- |
   | One sample | Supply `--bam`, `--vcf`, or prepared inputs directly |
   | Several samples | Use `--input_csv samples.csv` and follow the sample-sheet rules |

For example, an Illumina WES BAM uses `--bam sample.bam --type short`. Add
`--mode snp` for small variants only, add `--mode sv` for the structural and
repeat route, or leave `--mode` out of the command for combined analysis. There
is no `--mode combined` value.

For an existing VCF, use `--vcf` with the required `--mode snp` or `--mode sv`
and omit `--type`. For prepared ANNOVAR input and batch routing, use the
[input guide](docs/INPUTS.md) and [running guide](docs/USAGE.md).

## Quick start

### 1. Provide the required software and data

PipeVar documentation targets Linux execution through Docker or Singularity.
Install Nextflow with Java and a container runtime, obtain ANNOVAR through its
registration process, and prepare an hg38 reference for your data.
Aligned-read examples also require a discoverable BAM/CRAM index and
`<reference>.fai`; see [reference files](docs/INPUTS.md#reference-files) and
[alignment indexes](docs/INPUTS.md#alignment-indexes).

The configured first-attempt tasks request as much as 16 CPUs and 64 GB of RAM;
retries can request more. These are task allocations rather than a whole-run
minimum. See [environments and resources](docs/INSTALLATION.md#execution-environments)
and the [30× example benchmark](docs/BENCHMARKS.md).

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
