# Module conventions

The active tree contains 41 module leaves and 45 process definitions. A module
represents a tool operation, not execution cardinality; there are no active
`multi_*` module directories or process declarations.

## Contracts

- Sample-specific inputs and outputs carry `sample_id` as tuple element zero.
- References and run-wide settings use separate value channels.
- Reused outputs have named emits.
- Optional file and region inputs use a staged path or `[]`, never `"null"`,
  `"no"`, or another sentinel string.
- Modules do not infer sample identity from filenames.
- Containers remain beside tool commands. Resources, retries, publishing, GPU
  options, tool policy, and executor settings remain in `nextflow.config`.

For example, `(sample_id, BAM, BAI)` is a sample-keyed alignment record: every
downstream join matches element zero, while the remaining fields keep their
documented order. Internal `small_variant` names include SNVs and indels; the
public `snp` vocabulary remains available for CLI compatibility.

Ownership is explicit:

- 32 PipeVar tool/semantic leaves remain directly under `modules/`.
- `modules/local` contains the `legacy_artifact_compat` path for the
  `MATERIALIZE_PUBLIC_ARTIFACT` compatibility process, the validated shared
  `gatk_resource_bundle`, and PipeVar-specific `longphase_support`.
- `modules/nf-core` contains nine unmodified modules locked by `modules.json`:
  DeepVariant, ExpansionHunter, Manta, Sniffles, three GATK stages, and two
  LongPhase stages.

## Intentional multi-process modules

Distinct behavior or incompatible outputs justify separate processes even when
the executable is shared. `ngs_prio`, `snp_prio`, `sv_prio`, and
`longphase_support` each expose one age-capable prioritizer for both single and
CSV inputs; single inputs pass a blank age value. The report module contains
keyed no-repeat, repeat, and repeat-plus-mito operations. These are biological
or compatibility boundaries, not single-versus-multi duplicates.

The former local HaplotypeCaller/GATK-resource and LongPhase composites are
replaced by shared keyed helpers and unmodified upstream stages. CNVnator and
mitochondrial modules remain intact. Historical cardinality-specific module paths and process names are removed
internal APIs. Public CLI behavior, output filenames, and publishing contracts
remain supported, but renamed processes do not retain practical compatibility
with old `-resume` task state.
