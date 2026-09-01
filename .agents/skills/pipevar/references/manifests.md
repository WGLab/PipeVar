# Sample-sheet manifests

Use this reference whenever `--input_csv` is present. PipeVar supports a basic sheet for one shared input family and a typed sheet for explicit row-level paths. See [INPUTS.md](../../../../docs/INPUTS.md) for maintained examples.

## Choose one schema

| Schema | Use when | Route selector |
| --- | --- | --- |
| Basic | Every row uses the same BAM/CRAM or VCF family | Launch with `--bam true` or `--vcf true` |
| Typed | Rows explicitly declare supported input kinds and path fields | `input_kind` selects the family; omit legacy selector booleans |

Do not mix the two launch styles. Presence of an `input_kind` header selects the typed schema.

## Basic sheet

Required columns:

```csv
sample,file_path,note_path
```

Optional columns are `age_of_onset` and the legacy alternate `age`. The important non-obvious rule is that `note_path` can contain either clinical notes or HPO files:

- default behavior interprets `note_path` as a clinical note;
- launch with `--note no` when every `note_path` is an HPO file.

Examples of route-level selectors:

```text
Aligned sheet: --input_csv samples.csv --bam true --ref_fa hg38.fa --type short
VCF sheet:     --input_csv samples.csv --vcf true --ref_fa hg38.fa --mode snp
```

Replace type/mode for the actual route. Sample identifiers must be unique. Validate that every path exists and every row follows the shared family and phenotype interpretation.

## Typed sheet

Use this documented shared header for maximum interoperability across current advanced typed layouts:

```csv
sample,input_kind,phenotype_path,phenotype_format,age_of_onset,snv_txt_path,snv_vcf_path,sv_vcf_path,vcf_path,alignment_path,alignment_index_path
```

For ordinary documented runs and proband rows, required base fields are `sample`, `input_kind`, `phenotype_path`, and `phenotype_format`. Supported phenotype formats are `clinical_note` and `hpo`.

The under-documented `--denovo_filter yes` implementation permits father, mother, and sibling control rows to leave phenotype path/format blank; probands still require phenotype data. Do not rely on this exception unless de novo filtering was explicitly requested and its current code path has been inspected.

| `input_kind` | Required main paths | Run-level behavior |
| --- | --- | --- |
| `bam_ngs`, `cram_ngs` | `alignment_path`; index may be explicit or discoverable | Provide `--ref_fa` and explicit `--type`; optional mode |
| `vcf_snv`, `vcf_sv` | `vcf_path` | Provide `--ref_fa`; omit `--type` and `--mode` |
| `annotated_snv` | `snv_txt_path` and `snv_vcf_path` | Optional SV/alignment columns determine the supported prepared route |

The code requires the four base columns plus the path columns needed by the selected route; unrelated conditional columns may be absent. Prefer the shared header above when creating reusable manifests, leaving unused path cells blank.

## Row consistency

Before command generation, enforce these documented rules:

- `sample` values are nonempty and unique; prefer filename-safe values. Filename-safe characters are currently enforced when de novo filtering is enabled.
- Rows use a supported homogeneous input-kind group; BAM and CRAM aligned rows may be combined.
- Non-annotated typed sheets use one phenotype format per run.
- All `annotated_snv` rows consistently provide or omit an alignment.
- All `annotated_snv` rows consistently provide or omit `sv_vcf_path`.
- Alignment-backed annotated rows currently support short reads.
- Mitochondrial analysis requires an alignment for every relevant annotated row.
- VCF-only and annotation-only rows cannot produce read-backed repeat results.

For both basic and typed sheets, preflight uniqueness explicitly; do not assume every route validates duplicates before execution. Do not infer missing path fields from similarly named files unless the user asks to generate or repair the sheet. Ambiguous sample-to-file matches require review.

## Age and other clinical metadata

`age_of_onset` is optional. Accepted documented compact values are an integer or an integer followed by `d`, `m`, or `y`; bare integers are interpreted as years. Keep clinical metadata row-scoped.

Do not print clinical-note contents while validating a sheet. In summaries and bug reports, minimize or redact sample identifiers and absolute clinical-data paths.

## Generate a sheet

The repository helper is interactive:

```bash
./scripts/generate_input_csv.sh
```

It can create either schema and can add `sv_vcf_path` to a typed annotated-SNV sheet. Running it writes a new CSV and therefore requires user intent to create or update a manifest.

Review its output before launch because it:

- scans only the selected directory's top level;
- uses exact, normalized, containment, and first-token matching;
- may skip ambiguous or incomplete pairs;
- may leave a missing index field blank;
- does not support commas or quoted CSV fields; and
- keeps the first duplicate inferred prefix with a warning.

Its update action uses a separate output path and does not overwrite the source sheet.

## Typed-sheet command decision

Do not answer a generic “run this typed sheet” request with one definitive command until these are known:

1. the sheet's `input_kind` set;
2. phenotype format;
3. sequencing type if any row is alignment-backed;
4. requested SNP, SV, or combined behavior where the route permits a mode;
5. execution profile, reference, and output directory.

Read [route-selection.md](route-selection.md) after inspecting the header and representative non-sensitive row metadata.
