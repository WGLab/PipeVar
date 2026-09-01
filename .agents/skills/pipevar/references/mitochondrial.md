# Mitochondrial analysis

Load this reference only for mitochondrial requests or mitochondrial result review. PipeVar runs mtDNA analysis alongside supported alignment-driven nuclear routes and publishes mitochondrial candidates separately.

## Eligibility

Enable with:

```text
--mito yes
```

Required route conditions:

- input includes BAM/CRAM alignment data;
- analysis is SNP-only (`--mode snp`) or supported combined analysis (omit `--mode`);
- `--mode sv` is not present; and
- long-read input does not use `--light yes`.

VCF-only and annotation-only inputs without an alignment cannot run mitochondrial analysis. Typed prepared inputs require an alignment path for every row participating in mitochondrial analysis. Current direct CLI validation rejects `--mito yes` whenever `--vcf` is also set, so a single prepared-ANNOVAR command cannot add mtDNA even if it supplies `--bam`; use a supported typed annotated manifest with alignment instead.

## Short-read route

With `--type short`, PipeVar prepares mitochondrial reads and uses Mutect2 in mitochondrial mode. The reference bundle must include:

- `<reference>.fai`;
- `<reference-basename>.dict`; and
- `<reference>.amb`, `.ann`, `.bwt`, `.pac`, and `.sa` BWA sidecars.

Confirm the alignment and reference use compatible contig names and sequence content. Do not create missing sidecars unless the user asked to prepare the reference.

## Long-read route

With `--type ont` or `--type pacbio`, PipeVar uses a mitochondrial Clair3 route followed by postprocessing. It requires the reference FASTA index but not the short-read BWA sidecars or sequence dictionary.

The route uses the standard long-read caller, so `--light yes` is invalid because light mode selects NanoCaller.

## Mitochondrial contig

`--mito_contig` supplies the preferred alias and defaults to `chrM`. Alias fallback is not uniform across stages: BAM extraction recognizes `MT`, `chrM`, `M`, and `chrMT`, while the current short-read Mutect2 reference fallback checks only `MT` and `chrM`.

Before launch, compare:

- FASTA `.fai` contigs;
- BAM/CRAM index contigs; and
- any supplied mitochondrial alias override.

Require an explicit `--mito_contig` that matches both BAM/CRAM and FASTA, especially for `M` or `chrMT`; do not rely on end-to-end fallback for those aliases. Alias selection cannot fix a reference/alignment sequence mismatch. If multiple mitochondrial-like contigs exist, require an explicit decision rather than guessing.

## Annotation and prioritization

Both short- and long-read routes normalize and index the mitochondrial callset, annotate it with mitochondrial evidence sources, and create a prioritized TSV.

In the current module, `--hmtvar_data` is only a truthy/falsy status marker passed to annotation as `bundled` or `asset_missing`; the file/path is not staged or read. There is no documented content/path contract or canonical enabled value. Omit it by default. If explicitly set, explain that any nonempty value currently changes only the status marker, record that value, and do not claim HmtVar evidence was consumed.

Two threshold groups are applied sequentially while constructing `*.mito.prioritized.tsv`:

- initial call-quality filters: `--mito_min_vaf`, `--mito_min_depth`, and `--mito_min_alt_reads`; and
- secondary strict priority gates: `--mito_gui_min_af`, `--mito_gui_min_apogee2`, and `--mito_gui_min_mitotip`.

Current code requires the secondary AF, APOGEE2, and MitoTip comparisons to pass together before a row is retained. Despite the `gui` parameter names and current documentation wording, these are not display-only settings; the report consumes the resulting prioritized table. Do not substitute one threshold group for the other. Use [PARAMETERS.md](../../../../docs/PARAMETERS.md) for configured defaults, verify semantics against the checked-out prioritization code, and report non-default thresholds with the run command.

## Expected artifacts

Use [OUTPUTS.md](../../../../docs/OUTPUTS.md) for the exact contract. Typical mitochondrial artifacts include:

- normalized/indexed `*.mito.vcf.gz` and `.tbi`;
- short-read duplicate metrics;
- `*.mito.annotated.tsv`;
- annotated/indexed mitochondrial VCF; and
- `*.mito.prioritized.tsv` as the primary review table.

The mitochondrial prioritized table remains separate from nuclear `.prio.vcf` files. A supported combined HTML report may incorporate mitochondrial evidence, but do not infer that the nuclear VCF contains mtDNA candidates.

## Completion checks

1. Confirm the expected mitochondrial VCF and index exist and are nonempty.
2. Confirm annotation and prioritized TSV headers are present even when no candidates pass.
3. Check that the resolved contig matches the reference and input alignment.
4. Record both mitochondrial threshold groups and the HmtVar status marker; do not report HmtVar data as used.
5. Review warnings from preparation, caller, normalization, and annotation tasks.

Mitochondrial ranks and scores are evidence for expert review. Do not convert them into pathogenicity claims or diagnoses without independent clinical interpretation.
