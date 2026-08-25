# PipeVar output guide

PipeVar publishes selected files directly into `--output_directory`, which
defaults to the launch directory. Exact files depend on the input, sequencing
type, analysis mode, and enabled documented features.

## Start with these results

For a combined route, open `*.variant_html_report.html` first; it is the easiest
human-readable summary. For a route without an HTML report, start with
`*.prio.vcf`. For mitochondrial review, start with
`*.mito.prioritized.tsv`.

In the patterns below, `*` is the sample prefix. VCF is structured variant data,
while TSV files are tab-separated tables that can be opened in spreadsheet or
data-analysis software.

| Pattern | Meaning |
| --- | --- |
| `*.prio.vcf` | Prioritized candidate variants |
| `*.prio_gene.vcf` | Gene-focused prioritized candidates |
| `*.frequency_audit.tsv` | Frequency source, threshold, decision, and evidence audit |
| `*.variant_html_report.html` | Integrated report for supported combined/report routes |
| `*.rank_var.tsv` | RankVar small-variant evidence |
| `*.rankscore_filtered.tsv` | RankScore-filtered small-variant evidence |
| `*.clinvar.txt` | Accepted ClinVar pathogenic/likely-pathogenic evidence |
| `*.mito.prioritized.tsv` | Prioritized mitochondrial candidates when enabled |

These files contain candidates for expert review rather than a diagnosis.

## Small-variant calls and annotations

| Pattern | Condition | Meaning |
| --- | --- | --- |
| `*.deepvariant.vcf.gz` | Short-read default route | DeepVariant calls |
| `*.recal.vcf.gz` | Supported short-read light route | HaplotypeCaller/VQSR calls |
| `*.clair3.vcf.gz` | Long-read default route | Clair3 calls |
| `*.nanocaller.vcf.gz` | Supported long-read light route | NanoCaller calls |
| `*.hg38_multianno.txt` | ANNOVAR small-variant annotation | hg38 annotation table |
| `*.hg38_multianno.vcf` | ANNOVAR small-variant annotation | Annotated VCF where published |

Nuclear caller VCF indexes are not consistently published and are not part of
the documented output contract.

## Structural and copy-number outputs

| Pattern | Condition | Meaning |
| --- | --- | --- |
| `*_manta.vcf` | Short-read structural route | Manta calls |
| `*_cnvnator.tab`, `*_cnvnator.vcf` | CNVnator enabled | Copy-number calls |
| `*.shortread_sv.merged.vcf` | Multi-caller short-read route | Merged evidence before collapse |
| `*.shortread_sv.truvari_collapsed.vcf` | Multi-caller short-read route | Collapsed structural-variant evidence |
| `*.sniffles.vcf` | Long-read structural route | Sniffles calls used for analysis and phasing context |
| `*.exonic.vcf` | Structural annotation | ANNOVAR-filtered exonic records |
| `*.phenosv.filtered.tsv` | PhenoSV scoring | Passing phenotype-scored records |

When common structural-variant filtering is enabled, PipeVar also publishes:

| Pattern | Meaning |
| --- | --- |
| `*.common_sv_filtered.vcf` | Unmatched records and matches retained below the upstream ceiling |
| `*.common_sv_removed.vcf` | Matches removed above the upstream ceiling |
| `*.common_sv_filter.summary.tsv` | Match and retain/remove counts |

## Repeat-analysis outputs

| Pattern | Condition | Meaning |
| --- | --- | --- |
| `*.json` | Short-read aligned reads | ExpansionHunter result JSON |
| `*.eh.tsv` | Short-read aligned reads | Filtered ExpansionHunter table |
| `*_nanorepeat_result.tsv` | Long-read aligned reads | Final NanoRepeat comparison table |

VCF-only and annotation-only routes do not produce read-backed repeat results.

## Mitochondrial outputs

| Pattern | Meaning |
| --- | --- |
| `*.mito.vcf.gz`, `*.mito.vcf.gz.tbi` | Normalized and indexed mitochondrial calls |
| `*.mito.dup_metrics.txt` | Short-read mitochondrial duplicate metrics |
| `*.mito.annotated.tsv` | Mitochondrial annotation table |
| `*.mito.annotated.vcf.gz`, `*.mito.annotated.vcf.gz.tbi` | Annotated and indexed mitochondrial VCF |
| `*.mito.prioritized.tsv` | Candidates passing mitochondrial prioritization |

Mitochondrial candidates remain separate from nuclear `.prio.vcf` files and can
be included in the HTML report when the selected route creates one.

## Phenotype outputs

| Pattern | Condition |
| --- | --- |
| `*_phenotagger_patient_hpo.txt` | Clinical note processed by PhenoTagger |
| `*_phenogpt2_patient_hpo.txt` | Clinical note processed by PhenoGPT2 |
| `*_phen2gene*` | Phen2Gene rankings and optional targeted-region files |

## Published results and the work directory

The Nextflow `work/` directory contains staged inputs and additional
intermediates. Preserve it when `-resume` or task-level debugging is required.
For routine review and transfer, start with the published files in
`--output_directory`.

Return to [Running PipeVar](USAGE.md), the [workflow guide](WORKFLOW.md), or the
[documentation map](../README.md#documentation).
