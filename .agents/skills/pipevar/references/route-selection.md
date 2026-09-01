# Route selection

Use this reference to translate the user's starting material and analysis goal into a supported PipeVar route. The detailed operator examples remain in [USAGE.md](../../../../docs/USAGE.md), and exact parameters remain in [PARAMETERS.md](../../../../docs/PARAMETERS.md).

## First classify the input

| Input | Selection | Required route facts | Key restrictions |
| --- | --- | --- | --- |
| Single BAM/CRAM | `--bam <PATH>` | phenotype, `--ref_fa`, explicit `--type`, desired branch | Alignment index required |
| Single VCF | `--vcf <PATH>` | phenotype and explicit `--mode snp` or `--mode sv` | Omit `--type`; reference rules differ below |
| Prepared ANNOVAR SNV | `--annotated_snv yes --annovar_txt <TXT> --vcf <VCF>` | phenotype and matching prepared pair | SNP-led; advanced combinations use a typed sheet |
| Basic batch sheet | `--input_csv <CSV>` plus `--bam true` or `--vcf true` | one shared input family and phenotype interpretation | Legacy selector is mandatory |
| Typed batch sheet | `--input_csv <CSV>` | row `input_kind`, shared route facts | Do not add `--bam true` or `--vcf true` |

Do not treat FASTQ as a PipeVar input. PipeVar begins with aligned BAM/CRAM, an existing VCF, or prepared annotations.

## Select the analysis branch

| Goal | Aligned-read invocation | Consequence |
| --- | --- | --- |
| Small variants only | `--mode snp` | Small-variant calling/annotation/prioritization |
| Structural/copy-number variants and repeats | `--mode sv` | SV/CNV and read-backed repeat route |
| Combined nuclear analysis | Omit `--mode` | Supported small-variant and SV/CNV/repeat integration; supported routes add an HTML report |

Never write `--mode combined`. For a single VCF, choose exactly `--mode snp` or `--mode sv`; supplied calls cannot select the aligned-read combined route.

## Select sequencing type and callers

Always set a sequencing type for aligned reads.

| `--type` | Default small-variant caller | `--light yes` caller | Structural caller | Repeat caller |
| --- | --- | --- | --- | --- |
| `short` | DeepVariant | GATK HaplotypeCaller | Manta plus CNVnator when enabled | ExpansionHunter |
| `ont` | Clair3 | NanoCaller | Sniffles | NanoRepeat |
| `pacbio` | Clair3 | NanoCaller | Sniffles | NanoRepeat |

Use light mode with the documented small-variant-only route. Do not infer that a supplied VCF or prepared annotation can gain read-backed repeat or mitochondrial analysis without an alignment.

## VCF reference rules

- Single VCF: `--ref_fa` is recommended and becomes required for `--target yes`.
- Basic or typed VCF sample sheet: provide `--ref_fa`.
- Omit `--type` for VCF-only routes.
- Typed `vcf_snv` and `vcf_sv` rows infer the mode from `input_kind`; do not add a contradictory run-level `--mode`.

## Prepared-annotation routes

A prepared small-variant input consists of a paired hg38 multianno TXT and VCF. Supported typed-sheet layouts include:

- annotated SNV only;
- annotated SNV plus a supplied annotated SV VCF;
- annotated SNV plus an alignment for called SV/repeat analysis; or
- annotated SNV plus supplied SV evidence and an alignment.

Prepared-SNV routes are SNP-led: use `--mode snp` or omit `--mode`, never `--mode sv`. `--target yes` is unsupported. If an annotated route includes an alignment, it currently supports short reads only. Without an alignment, do not promise read-backed repeats or mitochondrial analysis.

The documented advanced batch combinations use typed sheets. Direct single-sample annotated-SNV plus BAM or supplied-SV routes exist in the code but are not fully covered by the operator guides; inspect their current validation and wiring before proposing them. Current top-level validation also rejects `--mito yes` whenever a direct command sets `--vcf`, even if it also supplies `--bam`. For prepared-input mtDNA, use a supported typed `annotated_snv` manifest with an alignment.

## Batch route rules

- A basic aligned sheet uses `--bam true`, an explicit `--type`, and optional `--mode` selection.
- A basic VCF sheet uses `--vcf true`, an explicit SNP/SV mode, and no `--type`.
- A typed aligned sheet contains `bam_ngs` or `cram_ngs` rows and uses one explicit run-level `--type`.
- A typed VCF sheet contains `vcf_snv` or `vcf_sv` rows and omits both `--type` and `--mode`.
- A typed annotated sheet follows the populated path columns and consistency rules in [manifests.md](manifests.md).

Do not build a definitive typed-sheet command until `input_kind`, phenotype format, host profile, and—when aligned—the sequencing type are known.

## Optional behavior compatibility

| Feature | Allowed | Reject or qualify |
| --- | --- | --- |
| `--mito yes` | BAM/CRAM with SNP-only or combined analysis | VCF-only, annotation-only without alignment, SV-only, long-read light |
| `--target yes` | Supported calling/annotation routes with reference | Prepared annotated-SNV routes; VCF without reference |
| `--light yes` | Documented SNP-only aligned routes | Long-read mitochondrial analysis; do not assume combined support |
| `--gene` | Final prioritization restriction | Validate symbols or a readable one-symbol-per-line file |
| PhenoGPT2 | Clinical-note input with configured GPU/models | HPO-only input does not need it; do not enable implicitly |

`--xtea yes` is visible in the pipeline help but is not fully covered by the operator guides. If explicitly requested, verify the current implementation constraints: BAM/CRAM input, `--type short`, and either `--mode sv` or omitted mode; reject SNP-only use. Typed VCF manifests can escape the current top-level validation even though their subworkflows do not run xTEA, so absence of a validation error is not proof of execution. Verify an aligned short-read subworkflow and the expected xTEA process/output in the trace. Do not advertise or enable it by default.

If a requested combination is not described here or in the operator documentation, inspect `main.nf` validation before proposing it and state that the route is not part of the documented interface.
