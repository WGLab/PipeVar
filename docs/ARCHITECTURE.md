# PipeVar_simplified architecture

PipeVar uses Nextflow DSL2. `main.nf` preserves the existing command-line
interface and selects a biological route. The named `INPUT_PREPARATION` helper
validates inputs and normalizes single-sample and CSV manifests into the same
sample-keyed records. Tool commands remain
in `modules/`, route orchestration remains in `subworkflows/`, and execution and
publication policy remains centralized in `nextflow.config`.

The active graph has 15 named workflows and 45 process definitions.
Structure is checked by `scripts/check_structure.sh`. Pinned upstream modules
and PipeVar-owned modules have separate ownership. When Nextflow is available,
`scripts/check_channel_contracts.sh` exercises reordered keyed joins, wide-record
flattening, batch binding reconstruction, and strict duplicate/mismatch failures.

## Channel contracts

Every sample-specific record carries the stable `sample_id` as tuple element
zero. Common shapes are:

| Record | Shape |
| --- | --- |
| Alignment | `(sample_id, alignment, index)` |
| Alignment with phenotype | `(sample_id, alignment, index, phenotype)` |
| VCF | `(sample_id, vcf)` |
| Annotated SNV | `(sample_id, annotation_txt, annotation_vcf, phenotype, phenotype_format)` |
| Clinical metadata | `(sample_id, age_of_onset)` |

References and run-wide controls remain separate value channels. Sample
associations use keyed joins with mismatch or duplicate checks; filenames must
not be used to reconstruct a dropped sample key. Optional regions have one
internal representation: a staged path when targeting is enabled and `[]` when
disabled. Strings such as `"no"` and `"null"` are not internal path values.

Single-sample, one-row CSV, and multi-row CSV inputs therefore traverse the
same route workflows. Prioritizers accept the same optional age field in every
route; single-sample inputs pass it blank. The input origin remains explicit
only for the batch-only LongPhase haplotag publication contract. Single-sample
values are unwrapped only at the final output boundary.

## Layers

1. `main.nf` owns route selection and passes prepared inputs to the biological
   workflows. `INPUT_PREPARATION` owns CLI compatibility, validation, reference
   preparation, and sample-record normalization.
2. `subworkflows/` owns biological route families and shared input-preparation,
   GATK, and LongPhase helpers.
   Cardinality-specific `SINGLE_*` and `INPUT_CSV_*` workflow APIs are removed.
3. `modules/` owns external tool operations. PipeVar operations remain in the
   existing flat directories, pinned upstream code is under `modules/nf-core`,
   and PipeVar-only adapters are under `modules/local`. Cardinality alone does
   not justify separate processes or directories.
4. `nextflow.config` owns resources, retries, publishing, runtime flags, bind
   mounts, GPU behavior, tool policies, and profiles.

## nf-core caller boundary

DeepVariant, ExpansionHunter, Manta germline, Sniffles, three GATK stages, and
two LongPhase stages are unmodified nf-core modules pinned by `modules.json`.
Routes translate
`(sample_id, ...)` records to nf-core `([id: sample_id], ...)` records only at
the caller boundary, then immediately restore the PipeVar sample-keyed shape.
Reference metadata is independent of sample metadata, and missing optional
files are always `[]`.

The local `MATERIALIZE_PUBLIC_ARTIFACT` process restores the public uncompressed
ExpansionHunter, Manta, and Sniffles filenames. Additional nf-core indexes,
gVCFs, candidate VCFs, realigned BAMs, SNF files, and version topics remain
internal. See [nf-core module provenance](NFCORE_MODULES.md) for the pinned
revision and migration boundary.

## Route ownership

The biological route families are `ngs`, `long`, `vcf`,
`preannotated_snp`, `preannotated_sv`, `preannotated_combined`, and
`mitochondrial`. `input_preparation`, `gatk_snp_calling`,
`longphase_processing`, and `preannotated_shared` are reusable helpers rather
than additional CLI routes. The `ngs` and `long` routes select SNP-only,
SV-only, or combined processing internally; `vcf` selects SNP or SV processing.
Full and light modes select callers inside the same active route;
preannotated supplied-SV and called-SV variants retain distinct wrappers but
reuse common evidence preparation; and mitochondrial long- and short-read
routes are colocated.

## Deliberately separate behavior

Three prioritization groups retain separate, purpose-named modules in a
shared module directory. Each group has one age-capable process for both
single-sample and CSV inputs; single-sample mode passes a blank age value.

GATK now acquires and validates one hg38 VQSR bundle per run and broadcasts it
to keyed HaplotypeCaller, VariantRecalibrator, and ApplyVQSR tasks. LongPhase
uses upstream preparation-independent phase and haplotag stages, while local
prioritization consumes the same optional age field in every route.
CNVnator and mitochondrial composites remain intentionally unchanged because
splitting them would add task boundaries without making their route contracts
clearer. Common-SV audit artifacts and report evidence combinations are also
preserved.

## Compatibility and resume boundary

The refactor preserves CLI parameters and aliases, defaults, manifest schemas,
analysis behavior, output filenames, publish locations, and top-level results.
Historical module paths, including the former local DeepVariant,
ExpansionHunter, Manta, Sniffles, HaplotypeCaller, GATK-resource-preparation,
and LongPhase paths, and the
`SINGLE_*`/`INPUT_CSV_*` named workflow APIs are internal APIs and are
intentionally not preserved.

Canonical process names and the rewired graph change task identities. Existing
published outputs remain compatible, but old `-resume` state from
`PipeVar_mito` or an earlier `PipeVar_simplified` graph is not practically
compatible with this graph; start a new run after upgrading.
