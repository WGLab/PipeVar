# Subworkflow conventions

The active source tree contains 15 named workflow definitions. Every route consumes sample-keyed records for both single and CSV
execution. `main.nf` performs input normalization and retains `batch_input`
only for the batch-only LongPhase haplotag publication contract.

| Directory | Owned workflows |
| --- | --- |
| `ngs` | mode-aware short-read SNP, SV, and combined route |
| `long` | mode-aware long-read SNP, SV, and combined LongPhase route |
| `vcf` | mode-aware supplied-VCF SNP and SV route |
| `preannotated_snp` | preannotated small-variant route |
| `preannotated_sv` | preannotated SNV plus supplied-SV route without alignment |
| `preannotated_combined` | supplied-SV and called-SV alignment-backed routes |
| `mitochondrial` | long-read and short-read mitochondrial routes |
| `input_preparation` | shared validation, normalization, and route selection |
| `gatk_snp_calling` | shared HaplotypeCaller/VQSR helper with one run-wide resource bundle |
| `longphase_processing` | shared phase, haplotag, and compatibility-sensitive prioritization helper |
| `preannotated_shared` | shared imported SNV/SV evidence preparation and strict combined-evidence assembly |

Sample-keyed means tuple element zero is the patient/sample ID; joins associate
records only when those IDs match and reject duplicates where the relationship
must be one-to-one. Avoid unrestricted `combine` on multi-sample channels
because it can pair every record from one channel with every record from the
other. Keep unrelated biological routes separate instead of hiding them
behind a universal optional-input tuple. Prioritizers receive an optional age
value in the same record contract for single and CSV inputs.

The historical `SINGLE_*` and `INPUT_CSV_*` workflow definitions, their include
paths, and the obsolete light-adapter names are removed internal APIs. Their CLI
functionality is preserved by the keyed routes, including `--light`, filenames,
publishing, and top-level results. Named internal emits are not a compatibility
contract. Rewired workflows and renamed processes are also a `-resume`
boundary; start a new run after upgrading.
