# Result review and completion checks

Use this reference to decide whether a run completed and which artifact to inspect first. Exact files are conditional on input, sequencing type, mode, and enabled features; [OUTPUTS.md](../../../../docs/OUTPUTS.md) is the maintained output contract.

## Choose the primary artifact

| Route | Start with | Caution |
| --- | --- | --- |
| Supported combined route | `*.variant_html_report.html` | Branch-only routes do not necessarily create HTML |
| SNP-only or SV-only route | `*.prio.vcf` | Caller-specific outputs are conditional for supplied VCF/annotations |
| Mitochondrial analysis | `*.mito.prioritized.tsv` | mtDNA candidates remain separate from the nuclear prioritized VCF |
| Gene-restricted review | `*.prio_gene.vcf` where produced | It is a restricted companion view, not a replacement for the full `.prio.vcf` |

Use sample prefix and route metadata to avoid selecting another sample's similarly named output.

## Route-aware supporting evidence

- Short-read default small-variant calling can publish a DeepVariant VCF; supported light routes use HaplotypeCaller/VQSR outputs.
- Long-read default small-variant calling can publish a Clair3 VCF; supported light routes use NanoCaller output.
- Short-read structural analysis can publish Manta, CNVnator, merged, and collapsed SV evidence depending on enabled callers.
- Long-read structural analysis uses Sniffles evidence.
- Short-read aligned routes can publish ExpansionHunter JSON/table results; long-read aligned routes can publish NanoRepeat results.
- VCF-only and annotation-only routes do not produce read-backed repeat results.
- Common-SV filtering, when enabled, can publish retained, removed, and summary files.
- `*.frequency_audit.tsv` records frequency source, threshold, decision, and evidence for final candidates where produced.

Do not list every possible artifact in a routine response. Report the primary result and supporting files relevant to the user's question.

## Minimal completion checks

1. Confirm Nextflow ended successfully and no earlier process failure is hidden by later log messages.
2. Confirm the expected primary artifact exists in `--output_directory` and is not truncated.
3. For VCFs, verify a valid header and expected sample identity; distinguish a valid header-only/zero-candidate result from a missing or malformed file.
4. For TSVs, verify the documented header exists even when no candidate rows pass.
5. For compressed/indexed mitochondrial VCFs, verify both the VCF and index are present.
6. Confirm report/trace outputs exist when enabled and retain them with the run provenance.
7. Review warnings about skipped, malformed, filtered, or unclassifiable records before declaring the biological analysis complete.

Useful read-only checks, when the corresponding tools exist, include:

```bash
test -s <PRIMARY_OUTPUT>
bcftools view -h <PRIORITIZED_VCF>
head -n 2 <PRIORITIZED_TSV>
```

Do not print patient-level records merely to prove that a file exists. Prefer headers, counts, and redacted summaries.

## Mitochondrial results

Expect a normalized/indexed mitochondrial callset, annotation outputs, and `*.mito.prioritized.tsv` when mitochondrial analysis completes. Short-read routes may also publish duplicate metrics. Apply the additional checks in [mitochondrial.md](mitochondrial.md).

Mitochondrial `mito_min_*` values are initial call-quality filters, while the currently named `mito_gui_*` values are secondary strict gates applied to the prioritized TSV before the report consumes it. Record both groups and read [mitochondrial.md](mitochondrial.md) for the code/documentation discrepancy.

## Interpret cautiously

PipeVar prioritizes candidate variants and genes for expert review. It does not establish a diagnosis or independently validate pathogenicity.

When summarizing results:

- separate pipeline completion from biological findings;
- report filters, thresholds, reference build, and evidence provenance;
- distinguish nuclear and mitochondrial results;
- avoid exposing patient identifiers or note content;
- state when zero candidates pass rather than calling the run failed; and
- recommend qualified clinical review for diagnostic interpretation.

For algorithmic meaning and evidence integration, read [WORKFLOW.md](../../../../docs/WORKFLOW.md). Do not infer scoring semantics from filenames alone.

## Preserve reproducibility

Retain or record:

- PipeVar revision and dirty state;
- complete command and profile;
- Nextflow/runtime versions;
- reference and database/model versions;
- report and trace files, which default to the Nextflow launch directory rather than necessarily `--output_directory`;
- non-default thresholds; and
- the `work/` directory while resume or audit is needed.

If an expected artifact is absent, switch to [troubleshooting.md](troubleshooting.md) before rerunning.
