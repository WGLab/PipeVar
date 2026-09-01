# Input and environment preflight

Run only the checks relevant to the selected route. These commands inspect state; they do not prepare indexes or modify inputs. Use [INPUTS.md](../../../../docs/INPUTS.md) for the complete input contract.

## Host tools and effective configuration

From the PipeVar repository root, check the required host commands:

```bash
command -v nextflow
java -version
command -v docker
command -v singularity
```

Only the selected container runtime must be present. If Nextflow is available, inspect the selected profile without launching the workflow:

```bash
nextflow -version
nextflow run main.nf --help
nextflow config -profile local_docker
```

Replace `local_docker` with the intended named profile. Do not assume `standard` is unchanged because `setup.sh` can rewrite it. Avoid pasting full effective configuration into public reports because it can expose host paths.

For Docker, `docker info` verifies daemon access. For Singularity, use `singularity --version`. On SLURM, also confirm scheduler access and that every absolute data/resource path is visible from compute nodes.

## General file checks

Use absolute, quoted paths. Confirm inputs are readable and the output parent exists and is writable:

```bash
test -r /absolute/path/to/input
test -d /absolute/path/to/output-parent
test -w /absolute/path/to/output-parent
```

Inspect an existing output directory before launch. Do not delete or reuse it automatically. Choose a new directory or continue only when its relationship to the intended run is understood.

## BAM and CRAM

- BAM accepts `<file>.bam.bai` or `<basename>.bai` where the route/helper discovers it.
- CRAM accepts `<file>.cram.crai` or `<basename>.crai`.
- CRAM decoding requires the exact matching reference.
- The pipeline expects aligned, indexed input. Do not claim sorting or reference compatibility from a filename alone.

When samtools is available, use read-only checks such as:

```bash
samtools quickcheck -v /absolute/path/sample.bam
samtools view -H /absolute/path/sample.bam
samtools idxstats /absolute/path/sample.bam
```

Use the header to confirm coordinate sort metadata, sample/read-group identity, and contig naming. Treat absent or contradictory metadata as an ambiguity to resolve, not permission to repair the file. `samtools quickcheck` is not a full biological or index validation.

## Reference FASTA

Alignment-driven routes require `--ref_fa` and `<reference>.fai`. Confirm both are readable:

```bash
test -r /absolute/path/hg38.fa
test -r /absolute/path/hg38.fa.fai
```

PipeVar's active ANNOVAR modules use hg38. Confirm the FASTA, alignments, supplied calls, indexes, and annotation databases describe the same build. Compare contig naming rather than assuming `hg38` and `GRCh38` resources are interchangeable in representation.

When samtools or bcftools is present, useful read-only views include:

```bash
cut -f1 /absolute/path/hg38.fa.fai
bcftools view -h /absolute/path/sample.vcf.gz | rg '^##contig='
```

Inspect representative contigs and mitochondrial aliases; do not emit entire sample-level variant content during preflight.

Short-read mitochondrial analysis additionally requires the sequence dictionary beside the FASTA and BWA sidecars:

```text
<reference-basename>.dict
<reference>.amb
<reference>.ann
<reference>.bwt
<reference>.pac
<reference>.sa
```

Long-read mitochondrial analysis requires the FASTA index but not those BWA sidecars or the sequence dictionary.

## VCF and prepared annotations

- Confirm the VCF is readable and structurally valid with an available VCF tool.
- Do not invent a `.tbi` requirement when the selected documented route does not state one.
- Confirm contigs/build agree with any supplied reference and hg38 annotation contract.
- Prepared ANNOVAR SNV input requires the paired multianno TXT and matching VCF; do not mix files from different annotation runs.
- Do not publish or copy full prepared inputs merely to validate them.

For a single VCF, a reference is recommended and required for targeted analysis. VCF sample-sheet routes require a reference.

## Phenotype inputs and privacy

Every analyzed sample needs exactly one phenotype source:

- `--hpo <FILE>` containing standard `HP:NNNNNNN` identifiers, normally one per line; or
- `--note <FILE>` containing a clinical note.

Direct clinical-note input uses PhenoTagger unless PhenoGPT2 is explicitly selected. Confirm files exist without printing note contents. When checking HPO formatting, report invalid line numbers or identifiers only; do not echo unrelated patient data.

Sample sheets carry phenotype paths and formats per the rules in [manifests.md](manifests.md). Treat sample identifiers, note paths, and logs as potentially identifying. Redact them in external issue reports unless disclosure is explicitly authorized.

## Pipeline resources

Inspect the effective `annovar_host_path` and `phenosv_host_path` in configuration. Both host directories must exist and be visible wherever tasks execute; profiles mount them at `/annovar` and `/PhenoSV/train_data`.

If a clinical-note run explicitly selects PhenoGPT2, apply the model, cache, GPU, and fresh-work-directory requirements in [runtime-and-resume.md](runtime-and-resume.md).

Do not run `setup.sh` as a preflight check. It downloads data and rewrites configuration.

## Preflight outcome

Before launch, summarize:

- resolved route and profile;
- input and phenotype paths checked;
- reference build and companion files checked;
- external resource mounts checked;
- output location and collision state;
- remaining warnings or assumptions.

Stop before execution if a required path, index, route choice, reference relationship, or authorization is unresolved.
