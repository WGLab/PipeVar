# PipeVar interpretability and reliability overhaul

## Goal

Make `PipeVar_simplified` reviewable by a bioinformatician and maintainable by a
Nextflow developer without changing public CLI parameters, supported input
schemas, biological decisions, filenames, publishing paths, or top-level
results. `PipeVar_mito` remains the untouched backup. Rewired tasks establish a
new practical `-resume` boundary.

## Findings that drive the work

- Route and module ownership is substantially clearer, but `main.nf` still
  mixes validation, channel construction, and route dispatch across roughly
  1,700 lines. Manifest parsing and the single/CSV dispatcher are now
  consolidated; separating validation from orchestration remains future work.
- Wide positional tuples obscure evidence identity, especially in combined and
  LongPhase workflows.
- Cardinality-based route comments and opaque variables have been replaced by
  sample-keyed and biological-purpose names; remaining long tuples still need
  named contracts.
- Most sample-key, file/index, frequency, and coordinate checks are necessary.
  Manifest validation and runtime channels now share one normalized
  row collection, and ambiguity checks fail early.
- Global retries are now limited to resource/termination failures, with a
  separate transient-download policy for the GATK resource bundle.
- Existing structure checks still overvalue exact counts and source spelling.
  The corrected GATK resource-argument defect demonstrates why command and
  channel-contract fixtures remain necessary.
- Critical prioritization algorithms run from container-root scripts whose
  source and version are not indexed in this repository.

## Implementation phases

### 1. Correctness and executable contracts

- Supply every VariantRecalibrator `--resource` flag together with its staged
  VCF filename and verify the four exact flag/file pairs.
- Replace global unconditional retry with a conditional retry for resource or
  termination exit statuses; give the shared resource downloader its own
  transient-failure retry policy.
- Keep strict keyed joins, required file/index checks, LongPhase coordinate
  validation, and resource checksum/size validation.
- Add focused fixtures that exercise generated caller arguments and named
  outputs. Stub success alone is not biological or command-line proof.

### 2. Input normalization and one dispatcher

- Parse each manifest once into normalized rows, then apply schema, phenotype,
  and clinical-metadata validation to those rows.
- Reject duplicate headers, header-only manifests, unsafe sample/output names,
  and missing legacy `file_path`/`note_path` fields at the input boundary.
- Centralize BAM/CRAM index resolution while preserving the current explicit
  index and `.bai`/`.crai` fallback order and messages.
- Normalize single and CSV inputs into the same keyed channel variables before
  route selection. Use one route-dispatch block; age-capable prioritizers receive
  a blank age for single inputs and normalized age for CSV rows.

### 3. Readable route contracts

- Document every workflow `take:` record shape and keep `sample_id` at tuple
  position zero.
- Use small named context maps inside orchestration where evidence tuples are
  wide, while keeping files visible as `path` inputs at process boundaries.
- Remove unused internal arguments such as SV overlap controls from SNP-only
  workflows after reference checks.
- Rename misleading variables and comments by biological purpose. Format long
  workflow calls vertically so reviewers can compare arguments.
- Remove the four unreachable `light.nf` definitions. Current `--light`
  behavior continues to select GATK or NanoCaller within the active routes.

### 4. Checks and documentation

- Treat module/workflow counts as reported inventory rather than the primary
  correctness gate. Continue checking ownership, empty directories, vendored
  integrity, resolvable includes, and forbidden cardinality-specific APIs.
- Replace fragile occurrence/text checks with single, one-row CSV, multi-row,
  missing-key, duplicate-key, targeted, age-metadata, and caller command
  fixtures.
- Add `workflow.md` as the short code-navigation guide. Keep
  `docs/WORKFLOW.md` for detailed biological behavior and correct its route
  claims about repeats and HTML reports.
- Add a provenance index for container-resident prioritization scripts,
  including input/output schemas and image versions.

## Acceptance criteria

- Nextflow 25.04.8 parses the complete DSL2 graph and all production profiles.
- Single, one-row CSV, and multi-row fixtures select the same keyed route for
  equivalent inputs; missing and duplicate sample keys fail before execution.
- GATK, DeepVariant, Manta, Sniffles, ExpansionHunter, and LongPhase command
  contracts retain expected arguments and filenames.
- Unified prioritizers accept age metadata for every route; single inputs pass a
  blank value and CSV rows retain optional age-of-onset. CNVnator, common-SV
  audit, mitochondrial, repeat, and HTML filenames are unchanged.
- `PipeVar_mito` matches its pre-work content checksum. No Git repository is
  created.

## Current implementation status

Completed in this overhaul:

- one normalized manifest row collection for validation and execution;
- one single/CSV dispatcher with shared age-capable prioritization;
- centralized alignment-index resolution and early input ambiguity checks;
- conditional retry policy and corrected GATK resource arguments;
- biological-purpose route variables and a named LongPhase evidence context;
- dynamic structure metrics, manifest boundary fixtures, and concise workflow
  navigation documentation.

Still required before claiming biological parity:

- representative-data comparisons for every caller and route;
- broader successful single/one-row/multi-row channel fixtures and mitochondrial
  branches;
- biological comparisons and image-digest capture for the container-resident
  scripts now indexed in `docs/CONTAINER_SCRIPTS.md`.

## Implementation order

1. Fix command correctness and retry behavior.
2. Add boundary validation and consolidate duplicated parsing helpers. (Done)
3. Unify route dispatch after input parity fixtures pass. (Done; verification ongoing)
4. Simplify names, tuple contexts, comments, and unused arguments. (In progress)
5. Update behavioral checks and documentation, then run the full verification
   matrix.
