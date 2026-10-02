# nf-core modules and container provenance

PipeVar vendors nine unmodified modules from `nf-core/modules` commit
`94b413835288b0e80bd29256d0a95e4fce020a8b`. `modules.json` is the component
lock even though this copied project intentionally has no Git metadata.

| PipeVar operation | Vendored process | Public PipeVar artifact |
| --- | --- | --- |
| DeepVariant | `DEEPVARIANT_RUNDEEPVARIANT` | `${sample}.deepvariant.vcf.gz` |
| ExpansionHunter | `EXPANSIONHUNTER` | `${sample}.json` |
| Manta germline | `MANTA_GERMLINE` | `${sample}_manta.vcf` |
| Sniffles | `SNIFFLES` | `${sample}.sniffles.vcf` |
| GATK HaplotypeCaller | `GATK4_HAPLOTYPECALLER` | internal `${sample}.vcf.gz` |
| GATK VariantRecalibrator | `GATK4_VARIANTRECALIBRATOR` | internal recalibration files |
| GATK ApplyVQSR | `GATK4_APPLYVQSR` | `${sample}.recal.vcf.gz` |
| LongPhase phase | `LONGPHASE_PHASE` | internal phased SNV/SV VCFs |
| LongPhase haplotag | `LONGPHASE_HAPLOTAG` | `${sample}_haplotag.bam` for CSV mode |

The vendored modules use nf-core metadata maps and standard named outputs.
PipeVar routes construct `[id: sample_id, single_end: false]` at each caller,
keep reference metadata separate, and unwrap `meta.id` immediately afterward.
The local `MATERIALIZE_PUBLIC_ARTIFACT` process decompresses caller outputs when
needed and materializes each artifact under PipeVar's required public filename.
Only those public compatibility artifacts are published; other
upstream outputs remain internal. `GATK_RESOURCE_BUNDLE` validates six objects
by MD5 and the composite dbSNP object by immutable size, then broadcasts the
bundle once per run. LongPhase VCF preparation and PipeVar-specific age-capable
prioritization remain local. The same prioritizer serves single-sample and CSV
inputs; single-sample records pass a blank age value.

DeepVariant remains on the public `--deepvariant_version` default of 1.9.0 and
selects the matching official CPU or `-gpu` image. ExpansionHunter retains its
content-addressed Wave image. Manta and Sniffles use the pinned Wave/BioContainers
images declared by their vendored modules. Dockerfiles in `docker_work` are not
rebuilt or rebased by this migration.

## Deferred candidates

Nuclear and mitochondrial Clair3 remain on the existing 1.2.0 implementation.
Clair3 2.0 must be evaluated separately for ONT/HiFi model mapping, targeted BED
behavior, genotype fields, filename compatibility, and downstream ranking
concordance; mitochondrial Clair3 is excluded from that evaluation.

Mitochondrial Mutect2, CNVnator, alignment normalization, Truvari, SURVIVOR,
xTea, and semantic PipeVar modules remain local. The mitochondrial and CNVnator
composites are intentionally not decomposed in this wave because the extra task
boundaries would not make their route contracts clearer.

NanoRepeat container work is also deferred. The active module expects
`/Nanorepeat_bed/nanorepeat.input.bed` and `compare_nanorepeat.sh`, while the
nominal `docker_work/nanorepeat` Dockerfile does not copy those files; they are
present only in the separate `nanorepeat_with_annotation` context. Reconcile
and contract-test that source before rebuilding its image.

## Updating vendored modules

Resolve a specific nf-core/modules commit, install all nine modules from that
same revision, and update `modules.json`. Do not edit installed files. Re-run
the structure, sensitive-contract, caller-boundary, configuration, and route
smoke checks before deleting the previous copies. Any module revision or caller
graph change establishes a new practical `-resume` boundary.
