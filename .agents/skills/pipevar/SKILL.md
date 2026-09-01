---
name: pipevar
description: "Operate PipeVar in this repository: choose supported alignment, VCF, prepared-annotation, batch, and mitochondrial routes; validate inputs and profiles; generate or execute commands; set up, resume, or troubleshoot runs; and review PipeVar output files and reports. Use when PipeVar, its CLI, sample sheets, profiles, logs, or reports are involved; do not use for generic Nextflow development or standalone clinical interpretation."
---

# Use PipeVar

Operate PipeVar from the repository root. Treat its checked-in documentation as the supported operator contract. Use `main.nf` and `nextflow.config` to confirm current implementation or mutable configuration, not to promote undocumented features.

## Operating workflow

1. Classify the request as explanation, command generation, execution, resume, troubleshooting, or result review.
2. Collect the run facts: input family, single or batch execution, phenotype format, sequencing type where applicable, analysis branch, runtime profile, reference, output directory, and optional mitochondrial, light, targeted, or GPU behavior.
3. If the route is unresolved, use [route selection](references/route-selection.md). Add [manifests](references/manifests.md) only for sample sheets and [mitochondrial analysis](references/mitochondrial.md) only for mtDNA work.
4. If a command must be built, use [command recipes](references/command-recipes.md). If the user supplied a route-complete command, validate it directly instead of loading recipes unnecessarily.
5. Before actual execution, apply [input preflight](references/input-preflight.md) and only the relevant sections of [runtime and resume](references/runtime-and-resume.md). Do not launch while a material choice or required input remains unresolved.
6. If the user requested guidance, return the validated command without launching it. If execution was explicitly requested, run after preflight without asking again unless new authorization or a consequential missing choice is required.
7. Load [result review](references/result-review.md) after completion. Load [troubleshooting](references/troubleshooting.md) only for a failure or resume decision.

## Load references by task

| Request | Required reference | Add when relevant |
| --- | --- | --- |
| Choose or explain a route | [route selection](references/route-selection.md) | [manifests](references/manifests.md), [mitochondrial analysis](references/mitochondrial.md) |
| Generate a command | [route selection](references/route-selection.md), [command recipes](references/command-recipes.md) | [runtime and resume](references/runtime-and-resume.md) for profile/GPU/resume choices |
| Execute a run | [input preflight](references/input-preflight.md), [runtime and resume](references/runtime-and-resume.md) | Route-specific references above |
| Create or validate a sample sheet | [manifests](references/manifests.md) | [input preflight](references/input-preflight.md) before launch |
| Diagnose or resume | [troubleshooting](references/troubleshooting.md) | [runtime and resume](references/runtime-and-resume.md) for cache decisions |
| Inspect results | [result review](references/result-review.md) | [mitochondrial analysis](references/mitochondrial.md) for mtDNA |

Do not load every reference automatically. Read only the route-specific material needed to answer or act.

## Universal route rules

- For aligned BAM/CRAM, always set `--type short`, `--type ont`, or `--type pacbio`; never rely on the configured `ont` default.
- Use `--mode snp`, `--mode sv`, or omit `--mode` for a supported combined aligned-read route. There is no `--mode combined`.
- A single VCF requires `--mode snp` or `--mode sv` and omits `--type`. A reference is recommended, required with `--target yes`, and required for VCF sample-sheet routes.
- Typed aligned rows require an explicit run-level `--type`. Typed `vcf_snv` and `vcf_sv` rows infer the branch from `input_kind`; omit both `--type` and `--mode`.
- A direct `--note` uses PhenoTagger unless the user explicitly selects and configures PhenoGPT2. Direct HPO input uses `--hpo`.
- Mitochondrial analysis uses `--mito yes`, requires a supported alignment route, rejects `--mode sv`, and rejects long-read `--light yes`. Current direct CLI validation rejects any route that also sets `--vcf`; prepared-input mtDNA therefore uses a supported typed annotated manifest with alignment.
- Nextflow options use one hyphen, such as `-profile` and `-resume`. PipeVar parameters use two, such as `--bam` and `--ref_fa`.

## Safety and authority

- Use absolute paths for data, references, resources, and published outputs. Confirm an existing nonempty output directory matches the user's intent before launch.
- Treat `setup.sh` as a download and configuration mutation. It rewrites resource paths and `profiles.standard`; run it only when setup/configuration work is requested.
- Prefer explicit named profiles in generated commands because `standard` can be rewritten locally.
- Use `-resume` only when code, inputs, parameters, containers, and mounted resource contents are unchanged. Preserve `work/` while debugging or resume matters.
- Do not repeatedly change parameters and relaunch. Diagnose the first failing process from `.nextflow.log`, trace/report data, and its `.command.*` files.
- Avoid reproducing clinical-note contents or unnecessary patient identifiers in commands, logs, summaries, and bug reports. Refer to files by path and redact excerpts.
- PipeVar produces candidates for expert review, not a diagnosis. Do not make clinical claims solely from ranks or scores.

## Handle missing facts

- Ask only for facts that change the route or safety of execution; derive discoverable paths and configuration from the repository or host.
- When the user has not supplied required paths, return a clearly labeled non-runnable template and enumerate every unresolved value.
- Distinguish a user choice from an implementation fact and from an assumption. Never silently substitute one input family, sequencing technology, phenotype format, profile, or analysis branch for another.
- If a requested route is unsupported, explain the exact incompatibility and the nearest documented alternative without launching it.

## Produce an operational handoff

For command generation or execution, state:

- the resolved input family, sequencing type, analysis branch, phenotype source, and profile;
- the complete runnable command, or a clearly labeled non-runnable template with every placeholder listed;
- prerequisites checked and any checks the current host could not perform;
- why `--mode` and `--type` are present or intentionally absent;
- expected primary artifacts and output directory; and
- any privacy, resource, or reproducibility caveat that changes the run decision.

For troubleshooting, identify the first failing process, evidence inspected, failure class, whether the run contract changed, and whether a fresh run or `-resume` is justified. Do not report a run as successful solely because files exist; verify the route-appropriate primary artifact and workflow status.

## Sources of truth

- Use [README.md](../../../README.md) for the project overview and documentation map.
- Use [INPUTS.md](../../../docs/INPUTS.md), [USAGE.md](../../../docs/USAGE.md), and [OUTPUTS.md](../../../docs/OUTPUTS.md) as the detailed operator contracts.
- Use [PARAMETERS.md](../../../docs/PARAMETERS.md) for public options and confirm mutable defaults in `nextflow.config`.
- Use [INSTALLATION.md](../../../docs/INSTALLATION.md) for setup and runtime details, [WORKFLOW.md](../../../docs/WORKFLOW.md) for analysis mechanisms, and [BENCHMARKS.md](../../../docs/BENCHMARKS.md) only for capacity planning.

Do not copy full parameter tables, output catalogs, module inventories, or external version claims into responses. Extract only what changes the requested decision.

The code may contain options outside the documented operator interface. When explicitly asked about one, inspect its validation and workflow wiring, label its documentation status, and verify prerequisites before proposing use. Do not advertise or enable it by default.
