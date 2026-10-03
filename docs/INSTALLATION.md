# PipeVar installation and setup

This guide covers software prerequisites, `setup.sh`, execution profiles,
container mounts, GPU configuration, runtime resources, troubleshooting, and
reproducibility. Start with the shorter setup path in the
[root README](../README.md) if you only need a local first run.

## On this page

- [Prerequisites](#prerequisites)
- [Execution environments](#execution-environments)
- [What setup does](#what-setupsh-does)
- [Setup examples](#setup-examples)
- [Profiles and bind paths](#execution-profiles)
- [GPU execution](#phenogpt2-and-gpu-execution)
- [Runtime resources](#runtime-resource-model)
- [Troubleshooting](#troubleshooting)
- [Reproducibility](#reproducibility)

## Prerequisites

Provide the following before running PipeVar:

- Nextflow 24.10.6 or newer and Java 17 or newer; see the
  [verified compatibility matrix](NEXTFLOW_COMPATIBILITY.md);
- Docker or Singularity;
- SLURM only when using a cluster profile;
- a licensed ANNOVAR installation containing `annotate_variation.pl`;
- storage for downloaded annotation and phenotype resources; and
- an hg38 reference FASTA with the companion files required by the selected
  workflow.

ANNOVAR must be obtained through its
[registration and download process](https://www.openbioinformatics.org/annovar/annovar_download_form.php).

## Execution environments

PipeVar documentation targets containerized Linux execution. The repository
provides these configured environments:

| Host environment | PipeVar profile | Repository status |
| --- | --- | --- |
| Linux workstation or server with Docker | `local_docker` | Configured and documented |
| Linux workstation or server with Singularity | `local_singularity` | Configured and documented |
| Linux SLURM cluster with Singularity | `standard` or `slurm_singularity` | Configured and documented |
| macOS, Windows, or WSL | None | Not documented or verified |

“Configured and documented” means the profile exists in `nextflow.config`; it
is not a claim that every distribution or runtime version has been tested. The
repository records a maintainer-verified Nextflow and Java matrix in
[Nextflow compatibility](NEXTFLOW_COMPATIBILITY.md), but does not certify every
Linux distribution, container runtime, scheduler, CPU architecture, or GPU.
Operators should record these versions with every run.

Containers package most analysis tools, but Nextflow, Java, the container
runtime, filesystem mounts, and any scheduler still operate on the host.
Run Nextflow from a writable directory and ensure `NXF_HOME` is writable. Use
`local_docker` or `local_singularity` on a workstation; the `standard` and
`slurm_singularity` profiles require a working SLURM installation.
PipeVar currently uses the v1 DSL2 syntax parser. Nextflow 26.04 defaults to
the newer v2 parser, so export `NXF_SYNTAX_PARSER=v1` before running that
release line.

## What `setup.sh` does

Run `setup.sh` from the PipeVar repository root. The script:

1. validates that the configured ANNOVAR directory contains
   `annotate_variation.pl`;
2. downloads the hg38 ANNOVAR databases used by PipeVar;
3. downloads and extracts full or light PhenoSV resources;
4. downloads the Phen2Gene knowledge-base assets;
5. adjusts paths inside the downloaded PhenoSV resources;
6. records ANNOVAR and PhenoSV host paths in `nextflow.config`; and
7. replaces `profiles.standard` with the selected execution backend.

The downloaded ANNOVAR set includes refGene, cytoBand, ExAC 0.3, avsnp147,
dbNSFP 4.7a, gnomAD 4.1 exome/genome, ClinVar 2024-09-17, and GTEx v8 eQTL/sQTL
resources.

The script does not install Nextflow, Java, Docker, Singularity, SLURM, or
ANNOVAR itself. It also does not install a reference FASTA/index bundle, custom
ExpansionHunter catalog, optional HmtVar data, or PhenoGPT2 models and caches.
It does pull the PhenoGPT2 and PhenoTagger containers for the selected runtime.
For Singularity profiles, it validates checksum-named SIFs, stores them in the
configured image directory, and writes their absolute paths to
`nextflow.config`.

## Setup examples

### Guided setup

```bash
bash setup.sh
```

The guided flow prompts for resource locations, bind sources, and the execution
backend.

### Local Docker

```bash
bash setup.sh --non-interactive \
  --profile=local_docker \
  --annovar-dir=/data/annovar \
  --phenosv-dir=/data/PhenoSV_model
```

### Local Singularity

```bash
bash setup.sh --non-interactive \
  --profile=local_singularity \
  --annovar-dir=/data/annovar \
  --phenosv-dir=/data/PhenoSV_model \
  --phenotype-image-dir=/data/pipevar-containers
```

### SLURM with Singularity

```bash
bash setup.sh --non-interactive \
  --profile=slurm_singularity \
  --annovar-dir=/shared/annovar \
  --phenosv-dir=/shared/PhenoSV_model \
  --phenotype-image-dir=/shared/pipevar-containers
```

### Light PhenoSV resources

`light` is a positional setup argument:

```bash
bash setup.sh light --non-interactive \
  --profile=local_singularity \
  --annovar-dir=/data/annovar \
  --phenosv-dir=/data/PhenoSV_model
```

This chooses the smaller PhenoSV resource bundle. Workflow light-mode behavior
is described in [Running PipeVar](USAGE.md#light-mode).

## Setup options

| Option | Meaning |
| --- | --- |
| `--non-interactive` | Use supplied values and defaults without prompts |
| `--profile=standard` | Configure the default SLURM/Singularity backend |
| `--profile=slurm_singularity` | Configure SLURM with Singularity explicitly |
| `--profile=local_singularity` | Configure local execution with Singularity |
| `--profile=local_docker` | Configure local execution with Docker |
| `--annovar-dir=<DIR>` | Existing ANNOVAR installation and database directory |
| `--phenosv-dir=<DIR>` | PhenoSV resource destination |
| `--annovar-bind=<DIR>` | Host bind source when it differs from the ANNOVAR install path |
| `--phenosv-bind=<DIR>` | Host bind source when it differs from the PhenoSV resource path |
| `--phenotype-image-dir=<DIR>` | Destination for validated PhenoGPT2 and PhenoTagger SIFs; defaults to `./containers` for Singularity |
| `--skip-phenotype-images` | Skip both phenotype-container pulls and leave existing SIF settings unchanged |

When bind options are omitted, their sources default to the corresponding
resource directories.

For `local_docker`, setup runs `docker pull` for both phenotype images and
clears any configured SIF paths. For Singularity profiles, configure
`SINGULARITY_TMPDIR` and `SINGULARITY_CACHEDIR` before setup when the runtime
defaults do not have enough local space. Existing validated `*.current.sif`
pointers in the image directory are reused, so routine setup reruns do not
download the large images again. Automatically staged filenames contain the
final SIF SHA-256; the manual digest-pinned procedure below remains available
when the source registry digest must also appear in the filename.

## Configuration changes and reruns

The setup script edits `nextflow.config`. Selecting a setup profile replaces
`profiles.standard`; it does not rewrite the separately named profile with the
same backend. For example, `--profile=local_docker` makes `standard` use local
Docker while leaving the named `local_docker` block intact.

Setup redownloads several resources on repeated runs. Treat rerunning it as a
configuration and download operation, not as a no-op. Keep a copy of intentional
local configuration changes before rerunning setup.

## Execution profiles

| Profile | Executor | Containers | Use case |
| --- | --- | --- | --- |
| `standard` | SLURM by repository default | Singularity | Default cluster configuration; setup can replace it |
| `slurm_singularity` | SLURM | Singularity | Explicit cluster execution |
| `local_singularity` | Local | Singularity | Workstation or single server |
| `local_docker` | Local | Docker | Docker-capable workstation or server |

Select a named profile with the Nextflow option `-profile`, for example
`-profile local_docker`. Complete runnable commands are provided in the
[root README](../README.md#3-run-one-sample-with-docker) and
[running guide](USAGE.md).

## Container bind paths

| PipeVar parameter | Container path | Purpose |
| --- | --- | --- |
| `--annovar_host_path` | `/annovar` | ANNOVAR installation and databases |
| `--phenosv_host_path` | `/PhenoSV/train_data` | PhenoSV model resources |

The host paths must exist wherever tasks execute. On a cluster, use identical
absolute paths on the submit and compute nodes. ANNOVAR paths are mounted only
for ANNOVAR tasks, and the PhenoSV path is mounted only for PhenoSV tasks;
unrelated routes do not require either resource directory.

## PhenoGPT2 and GPU execution

PhenoGPT2 is used for clinical notes only. A supported run requires:

```text
--phenotype_extractor phenogpt2
--GPU yes
--phenogpt2_model_host_path /absolute/versioned/new_model
```

With negation processing enabled, also provide:

```text
--phenogpt2_negation yes
--phenogpt2_negation_model_host_path /absolute/versioned/negation-model
--phenogpt2_embedding_model_host_path /absolute/versioned/embedding-model
```

`--phenogpt2_cache_host_path` optionally mounts a pre-created writable cache.
Model directories are mounted read-only. The configured backend adds `--nv` for
Singularity or `--gpus 1` for Docker. Scheduler GPU options come from
`--gpu_cluster_options`, whose default is `--gres=gpu:1`.

Keep model directories immutable and versioned. All paths must be available at
the same absolute location on every node that can execute the task.

## Pre-staging large phenotype-extraction images

PhenoGPT2 and PhenoTagger are large Docker images. Converting either image to
SIF during workflow launch can exceed the default pull timeout or temporary
filesystem capacity. The Singularity profiles accept pre-staged images without
changing the source Docker images:

```text
--phenogpt2_sif /shared/pipevar/images/phenogpt2_0.4_sha256-<digest>.sif
--phenotagger_sif /shared/pipevar/images/phenotagger_sha256-<digest>.sif
```

Both parameters are optional. Without them, PipeVar retains the module Docker
tags and lets Nextflow pull and convert the images. A supplied path is always
validated, even if that extractor is not selected: it must be absolute and
resolve to a readable regular file. PipeVar canonicalizes the path and passes a
properly encoded `file:///` URI to Singularity. Docker profiles reject SIF
parameters rather than silently pulling or ignoring a different image.

PipeVar allows two hours for remote Singularity pulls. Increasing the timeout
does not correct insufficient space, permissions, registry access, or a cache
that is not visible to compute nodes.

### Cache and temporary storage

Use a per-user shared location for Nextflow's reusable SIF cache and large
local scratch on the conversion host for Singularity's OCI cache and temporary
workspace:

```bash
export PIPEVAR_SHARED=/shared/users/<user>/pipevar-containers
export PIPEVAR_SCRATCH=/local-scratch/<user>/pipevar-singularity
export NXF_SINGULARITY_CACHEDIR="${PIPEVAR_SHARED}/nextflow-cache"
export SINGULARITY_CACHEDIR="${PIPEVAR_SCRATCH}/oci-cache"
export SINGULARITY_TMPDIR="${PIPEVAR_SCRATCH}/tmp"

mkdir -p "${NXF_SINGULARITY_CACHEDIR}" \
  "${SINGULARITY_CACHEDIR}" "${SINGULARITY_TMPDIR}" \
  "${PIPEVAR_SHARED}/images"
chmod 700 "${NXF_SINGULARITY_CACHEDIR}" \
  "${SINGULARITY_CACHEDIR}" "${SINGULARITY_TMPDIR}"
```

`NXF_SINGULARITY_CACHEDIR` must be visible at the same absolute path from the
Nextflow launch host and every SLURM compute node. Do not share a writable
runtime cache between users, and use filesystems that support atomic rename.
Check both free bytes and inodes before conversion. Reserve at least three times
the largest source-image size plus headroom in `SINGULARITY_TMPDIR`; the current
approximately 15 GB PhenoGPT2 image therefore needs at least 45 GB plus
headroom. If the site exposes Apptainer as `singularity`, use the corresponding
`APPTAINER_CACHEDIR` and `APPTAINER_TMPDIR` variables required by that runtime.

### One-time SIF staging

Resolve each mutable tag to a registry digest before staging. On a host with
Docker access, record the published digests:

```bash
docker image inspect --format '{{join .RepoDigests "\n"}}' \
  beoungl/docker_test:phenogpt2_0.4
docker image inspect --format '{{join .RepoDigests "\n"}}' \
  beoungl/docker_test:phenotagger
```

Pull each digest-qualified image to a unique temporary filename in the final
shared image directory. Replace the placeholders with the recorded digests:

```bash
singularity pull \
  "${PIPEVAR_SHARED}/images/.phenogpt2_0.4_sha256-<digest>.sif.partial" \
  docker://beoungl/docker_test@sha256:<digest>
singularity pull \
  "${PIPEVAR_SHARED}/images/.phenotagger_sha256-<digest>.sif.partial" \
  docker://beoungl/docker_test@sha256:<digest>
```

Validate each temporary image before publication:

```bash
singularity inspect "${PIPEVAR_SHARED}/images/.phenogpt2_0.4_sha256-<digest>.sif.partial"
singularity sif list "${PIPEVAR_SHARED}/images/.phenogpt2_0.4_sha256-<digest>.sif.partial"
sha256sum "${PIPEVAR_SHARED}/images/.phenogpt2_0.4_sha256-<digest>.sif.partial"
```

Repeat the checks for PhenoTagger. Confirm the final filename is unused, then
atomically rename the validated file on the same filesystem and make it
read-only:

```bash
mv "${PIPEVAR_SHARED}/images/.phenogpt2_0.4_sha256-<digest>.sif.partial" \
  "${PIPEVAR_SHARED}/images/phenogpt2_0.4_sha256-<digest>.sif"
chmod 0444 "${PIPEVAR_SHARED}/images/phenogpt2_0.4_sha256-<digest>.sif"
```

Stage PhenoTagger the same way. Never overwrite an existing SIF path: retain
the source registry digest, final SIF SHA-256, and `singularity inspect`
metadata, and use a new versioned path whenever contents change. This keeps
Nextflow `-resume` behavior auditable. Preserve `NXF_SINGULARITY_CACHEDIR` for
other containers that still use normal Nextflow caching.

## Runtime resource model

PipeVar configures up to three retries. Many CPU and memory requests scale with
the retry attempt, while time limits remain fixed. Local profiles retain these
requests and do not automatically scale them down for a workstation.

These are configured per-task allocations, not measured consumption or a
whole-run minimum. Local capacity must cover the largest tasks that Nextflow
runs concurrently. Disk requirements are not currently quantified; allow space
for references, downloaded databases, container images, published results, and
the Nextflow work directory.

Representative first-attempt allocations are:

| Workload | CPU | Memory | Time |
| --- | ---: | ---: | ---: |
| DeepVariant | 16 | 64 GB | 72 h |
| PhenoGPT2 | 16 by default | 64 GB | 12 h |
| Clair3 | 8 | 32 GB | 72 h |
| NanoCaller | 8 | 32 GB | 48 h |
| Manta | 4 | 32 GB | 48 h |
| CNVnator | 2 | 32 GB | 72 h |
| LongPhase | 4 | 32 GB | 24 h |
| ANNOVAR or repeat analysis | 1-2 | 16 GB | 12-72 h |
| Mitochondrial processes | 1-4 | 4-8 GB | 2-8 h |
| Final filters and reports | 1 | 4-8 GB | 1-8 h |

Use `--deepvariant_max_forks` to cap concurrent DeepVariant tasks on shared or
local systems. PhenoGPT2 defaults to one concurrent task.

Observed values from a supplied 30× example are kept separate from these
scheduler requests. See [Resources and 30× example benchmark](BENCHMARKS.md) for
the measurements, their limitations, and the metadata still needed for a
reproducible comparison.

## Troubleshooting

### Nextflow fails before parsing the workflow

Confirm that `java -version` reports Java 17 or newer and that Nextflow can
write to both the launch directory and `NXF_HOME`. An
`UnsupportedClassVersionError` reporting class-file version 61 with a runtime
maximum of 55 means Java 11 attempted to launch a Java-17 Nextflow build.

### Local runs try to submit to SLURM

The repository's `standard` profile uses SLURM by default. Select
`-profile local_docker` or `-profile local_singularity`, or intentionally replace
`standard` through setup.

### Containers cannot see ANNOVAR or PhenoSV

Check the configured host paths and confirm that they exist on every execution
node. The profiles mount them at `/annovar` and `/PhenoSV/train_data`.

### A local task requests too many resources

Reduce task concurrency, set `--deepvariant_max_forks`, or use a scheduler with
appropriate capacity. A retry can request more CPU and memory than the first
attempt.

### A resumed PhenoGPT2 run uses changed models

Use a fresh work directory whenever external model contents change. Mounted
model contents are outside normal Nextflow input staging.

## Reproducibility

- The four vendored nf-core caller modules are locked in `modules.json`; do not
  edit files below `modules/nf-core`. PipeVar-specific output adaptation lives
  below `modules/local`.
- Record the PipeVar revision, command, profile, parameter values, reference and
  database versions, and container identifiers.
- Keep reference, ANNOVAR, PhenoSV, and PhenoGPT2 resources in immutable,
  versioned locations.
- Use `-resume` only when code, inputs, mounted resources, and containers are
  unchanged.
- Preserve the Nextflow report and trace, which are enabled by default.

The exact upstream revision, caller containers, retained output contracts, and
deferred Clair3 evaluation are documented in
[nf-core module provenance](NFCORE_MODULES.md).

Return to the [documentation map](../README.md#documentation) or continue with
the [input guide](INPUTS.md).
