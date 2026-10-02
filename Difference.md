# PipeVar_mito and PipeVar_simplified comparison

## Scope and definition

This report compares the current filesystem snapshots of `PipeVar_mito` and
`PipeVar_simplified`. It does not assume that the directories represent two
formal releases from one Git history. `PipeVar_mito` was clean at commit
`a2215de8a55d718dc989291fc674f73abfe9c664`; the simplified directory was
reviewed as its current working tree.

Here, a **fundamental change** means a change to a supported route, biological
algorithm, caller parameter, threshold, published result, or interpretation of
the result. Process names, channel layouts, validation timing, and file
organization are described separately unless they alter observable behavior.

## Executive verdict

`PipeVar_simplified` is a substantial workflow and interface refactor. It is
not merely a reformatted copy. Its active graph, process identities, channel
contracts, validation, resource handling, and some container packaging differ.

The core supported biological routes, caller choices, mitochondrial algorithms,
principal filtering thresholds, repeat catalogs, and most effective container
versions are preserved. CSV batching and optional age-aware prioritization remain
supported across the shared route graph.

There are two confirmed differences that may matter to valid runs:

1. PacBio LongPhase now receives the explicit `--pb` platform argument. The old
   configuration supplied no LongPhase platform argument for PacBio.
2. A single-sample LongPhase haplotag BAM is no longer published by the
   simplified configuration; CSV batch haplotag BAMs remain published.

Manta, Sniffles, and LongPhase also moved to different container builds and
nf-core command wrappers while retaining their declared upstream tool versions.
These substitutions require matched output-concordance runs before biological
equivalence can be claimed.

## Architecture and workflow graph

| Area | PipeVar_mito | PipeVar_simplified | Assessment |
| --- | --- | --- | --- |
| Entry workflow | Monolithic `main.nf`, 1,724 lines | `main.nf`, 274 lines, plus the named `INPUT_PREPARATION` helper | Structural refactor |
| Single and CSV execution | Separate `SINGLE_*`, `INPUT_CSV_*`, and `multi_*` implementations | Shared sample-keyed route workflows | Structural refactor; reduces duplicated logic |
| Light mode | Separate light workflow files in addition to caller selection | Caller selection inside each active route | Structural refactor; caller choices are preserved |
| Channel identity | Mixed scalar, single, and batch layouts | Records consistently keyed by `sample_id` | Internal interface change; missing or duplicate joins now fail explicitly |
| Upstream modules | Primarily local wrappers | Nine pinned nf-core modules plus compatibility adapters | Process and packaging change |
| GATK resources | Downloaded inside single calling or staged per batch path | One validated bundle broadcast to all keyed samples | Operational change; HaplotypeCaller, VQSR resources, annotations, SNP mode, and 99.0 threshold are preserved |
| Resume state | Old process names and graph | New canonical process names and graph | Existing cross-version `-resume` state is incompatible |

Evidence: `PipeVar_mito/main.nf:1107-1135,1140-1713`,
`PipeVar_simplified/main.nf:80-263`,
`PipeVar_simplified/subworkflows/input_preparation/main.nf:1223-1284`, and
`PipeVar_simplified/docs/ARCHITECTURE.md:38-67,81-112`.

## Biological route comparison

| Route or evidence stage | Comparison | Verdict |
| --- | --- | --- |
| Short-read small variants | DeepVariant by default; GATK HaplotypeCaller with `--light yes` | Preserved |
| Long-read small variants | Clair3 by default; NanoCaller with `--light yes` | Preserved |
| Short-read SV/CNV | Manta, default CNVnator, optional xTea, and Truvari consolidation | Preserved route and parameters; Manta packaging changed |
| Long-read SV | Sniffles with supporting read names retained | Preserved route and arguments; Sniffles packaging changed |
| Short-read repeats | ExpansionHunter followed by the existing filter | Preserved; exact effective ExpansionHunter image retained |
| Long-read repeats | NanoRepeat on SV and combined routes | Preserved |
| Annotation and ranking | ANNOVAR, RankScore, RankVar, PhenoSV, inheritance-aware prioritization, and report evidence | Preserved for equivalent valid inputs |
| Age-aware prioritization | Single and CSV routes previously used separate process contracts | Unified age-capable prioritizers; single inputs pass a blank age and CSV rows retain optional age-of-onset |
| Prepared SNV/SV input | Reuses annotations and integrates supplied or newly called SV evidence | Preserved; route validation is stricter |
| Short-read mitochondrial | Mitochondrial preparation, Mutect2, normalization, annotation, and prioritization | Preserved |
| Long-read mitochondrial | Haploid-sensitive Clair3, postprocessing, annotation, and prioritization | Preserved |

The default nuclear thresholds remain the same: gnomAD AD `0.001`, gnomAD AR
`0.01`, RankScore `0.50`, RankVar `0.05`, PhenoSV `0.50`, GQ `20`, AD `15`, and
NanoCaller missing-GQ depth `20`. The mitochondrial prioritization defaults
remain VAF `0.01`, depth `50`, and alternate reads `5`.

Both repeat catalogs are byte-identical between the directories:

| Catalog | SHA-256 |
| --- | --- |
| `data/variant_catalog.json` | `40530d576842e8202a07576d94d75df35dc78b8199ae16b764883d55ddde602b` |
| `data/variant_catalog_grch38.json` | `e9019b3419e12261603c6517dfc7ee2563b93a8b60148b7f412846e9533af2e4` |

Mitochondrial process command bodies, output names, thresholds, and effective
images are retained. The simplified modules add `sample_id` to records and
replace duplicate single/batch orchestration; they do not introduce a different
mitochondrial caller or prioritization algorithm.

## Confirmed observable differences

### PacBio LongPhase platform selection

`PipeVar_mito/nextflow.config:230-233` sets an empty argument for LongPhase in
PacBio mode. `PipeVar_simplified/nextflow.config:188-201` passes `--pb` to
`LONGPHASE_PHASE`. This is a caller-parameter change and can alter PacBio
phasing. It should be treated as potentially biological until a matched PacBio
comparison confirms the expected effect.

### Haplotag publication

The old `longphase|multi_longphase` process published all declared outputs for
both single and batch runs (`PipeVar_mito/nextflow.config:488-497`). The
simplified `LONGPHASE_HAPLOTAG` selector publishes `*_haplotag.bam` only when
`meta.batch_input`
(`PipeVar_simplified/nextflow.config:471-479`). Prioritized VCFs and frequency
audits remain published. This is a user-visible output-availability change for
single-sample combined long-read runs.

### Validation and failure behavior

The simplified workflow rejects contradictory unified VCF `input_kind` and
`--mode` values, requires references consistently for raw VCF routes, rejects
fractional positive-integer options, validates mitochondrial/alignment
combinations, and uses strict sample-keyed joins. These changes affect invalid,
ambiguous, or incomplete commands by failing before task construction. They are
correctness changes rather than new biological analyses.

`PipeVar_mito` can let an explicit mode override the VCF kind in a unified
manifest (`main.nf:1611-1617`). The simplified input preparation rejects this
conflict (`subworkflows/input_preparation/main.nf:337-344`).

### Retry and resource-mount behavior

`PipeVar_mito` retries all process failures up to three times. The simplified
configuration retries only exit statuses `137`, `140`, and `143`, with two
retries (`nextflow.config:210-214`). Deterministic input or tool errors therefore
stop earlier.

The old profiles mounted ANNOVAR and PhenoSV resources globally and did not
quote complete mount arguments (`PipeVar_mito/nextflow.config:633-670`). The
simplified configuration validates canonical paths, quotes mounts, and attaches
them only to ANNOVAR, ANNOVAR-SV, and PhenoSV tasks. This changes runtime
isolation and error timing, not the algorithms used by those tools.

## Effective container comparison

The effective image is the image after `nextflow.config` overrides. A newer
image written in a vendored module is not effective when the configuration
overrides it.

| Stage | PipeVar_mito | PipeVar_simplified effective image | Assessment |
| --- | --- | --- | --- |
| DeepVariant | `google/deepvariant:1.9.0` or `1.9.0-gpu` | Same configurable default | Effective version preserved; vendored module's 1.10.0 declaration is overridden |
| GATK HaplotypeCaller/VQSR | `broadinstitute/gatk:4.5.0.0` | Same for HaplotypeCaller, VariantRecalibrator, and ApplyVQSR | Version and VQSR policy preserved; stages decomposed |
| ExpansionHunter | `community.wave.seqera.io/library/expansionhunter:5.0.0--389ada7e191a4fba` | Same | Exact effective Docker reference preserved |
| Sniffles | `beoungl/docker_test:sniffles_0.1` | `community.wave.seqera.io/library/sniffles:2.8.0--c25a97c10afa095a` | Build and wrapper changed; both available sources declare Sniffles 2.8.0 |
| Manta | `beoungl/docker_test:manta` | `community.wave.seqera.io/library/manta_python:0eb71149179b3920` | Build and wrapper changed; both available sources declare Manta 1.6.0 |
| LongPhase phase/haplotag | Combined `beoungl/docker_test:longphase_0.4.0` process | `quay.io/biocontainers/longphase:1.7.3--hf5e1c6e_0` | Stage split and packaging changed; both available sources declare LongPhase 1.7.3 |
| LongPhase prioritization | `beoungl/docker_test:longphase_0.4.0` | Same custom image for local prioritization stages | Prioritization script image retained |
| Mito Mutect2 preparation/calling | `beoungl/docker_test:mito_mutect2_0.1` | Same | Preserved |
| Mito Clair3 | `hkubal/clair3:v1.2.0` | Same | Preserved |
| Mito postprocessing | `beoungl/docker_test:mito_clair3_postprocess_0.1` | Same | Preserved |
| Mito annotation/prioritization | `beoungl/docker_test:mito_annotation_0.4.2` | Same | Preserved |
| Other equivalent active custom stages | Tag used by each old module | Same tag in simplified | No image-reference change found |

The unchanged group includes ANNOVAR, ANNOVAR-SV/Truvari helper, CNVnator,
common-SV filtering, ExpansionHunter filtering, NanoCaller,
NanoRepeat, Phen2Gene, PhenoGPT2, PhenoSV, PhenoTagger, RankScore, RankVar,
SURVIVOR, report generation, and xTea.

The local `/home/beoungle/docker_work` tree contains plausible Dockerfiles for
the custom images, including Manta, Sniffles, and LongPhase. It has no usable
commit history or image-build manifest. The mutable `beoungl/docker_test:*`
tags are not digest-pinned, so the local sources cannot prove what was used to
build the published images. Equal tags demonstrate equal references, not
immutable image identity.

## Removed files and inactive tools

`PipeVar_mito` contains local modules for cuteSV, CNVpytor, Scramble, generic
Truvari, and older merge helpers that are absent from the simplified tree.
Searches of `PipeVar_mito/main.nf` and its active subworkflows found no includes
or calls for cuteSV, CNVpytor, Scramble, or the old long/short merge modules.
Their removal is dead-code cleanup, not removal of a supported active route.
The active short-read Truvari merge remains present through
`truvari_shortread_sv_merge`.

## Difference classification

| Classification | Differences | Fundamental? |
| --- | --- | --- |
| Structural only | Shared single/CSV workflows, input-preparation workflow, sample-keyed tuples, removal of duplicated light workflows, GATK decomposition, compatibility adapters | No intended biological change; cross-version `-resume` is incompatible |
| Runtime and error handling | Strict validation, fail-on-mismatch joins, selective retry policy, canonical process-specific mounts | Observable failure behavior; not a biological algorithm change for valid inputs |
| User-visible output | Single-sample LongPhase haplotag BAM is no longer published | Yes, for output availability |
| Potentially biological | PacBio LongPhase gains `--pb`; Manta, Sniffles, and LongPhase use different builds/wrappers | Requires matched concordance testing |
| Preserved biological behavior | Route matrix, caller choices, thresholds, mitochondrial commands, repeat catalogs, GATK VQSR policy, annotation and prioritization stages | No fundamental change found by static inspection |
| Removed inactive code | cuteSV, CNVpytor, Scramble, and unused merge wrappers | No active-route effect found |

## Limits of this comparison

Static inspection establishes intended parity and confirms the differences
above, but it cannot establish byte-identical or biologically identical output.
That requires matched fixtures for:

- short-read SNP, SV, and combined analysis;
- ONT and PacBio SNP, SV, and combined analysis;
- short- and long-read mitochondrial analysis;
- raw VCF and prepared-annotation routes;
- single-sample and multi-row CSV manifests; and
- CPU/GPU modes where applicable.

Those runs should compare called records, filters, INFO/FORMAT fields,
prioritized variants and genes, frequency audits, report evidence, and published
filenames. Image digests should also be recorded before using output concordance
to make a release-equivalence claim.
