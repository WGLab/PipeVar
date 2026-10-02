# PipeVar resources and 30× example benchmark

This page separates PipeVar's configured task allocations from observed
measurements. Use the [installation guide](INSTALLATION.md#runtime-resource-model)
for current CPU, memory, time, retry, and concurrency settings.

## How to use these numbers

The tables below reproduce measurements supplied for a 30× example. Here,
“30×” is the supplied coverage label and ordinarily refers to roughly 30 reads
covering each assayed base on average; the assay and calculation method were
not supplied. This example is not a
minimum hardware specification, a performance guarantee, or one end-to-end run:
the table includes short-read and long-read callers that are alternatives in
different PipeVar routes. Do not add the rows to estimate a whole run. Nextflow
may also execute independent tasks concurrently.

The following benchmark metadata was not supplied and must be recorded before
the example is reproducible or can define a tested environment:

- whether 30× refers to WES, WGS, or another assay, plus input size and read
  technology;
- PipeVar revision, full or light route, reference, sample count, and enabled
  optional analyses;
- CPU model and core count, host RAM, GPU model and driver where applicable,
  and storage type;
- host OS, Nextflow, Java, container runtime, scheduler, and container versions;
- whether CPU is peak or average percentage, and what 100 represents;
- whether RAM is peak resident memory, average memory, or another measure;
- whether time is per-process elapsed wall time; and
- task concurrency and cache state.

CPU values greater than 100 may represent aggregated multi-core utilization,
but the measurement method has not been confirmed. RAM strings and rounded
times are retained as supplied rather than converted under an assumed unit
definition.

## Reported measurements

### Shared annotation and prioritization

| Software | Role | CPU usage (reported) | RAM usage (reported) | Time (reported) |
| --- | --- | ---: | ---: | ---: |
| Phen2Gene | Phenotype ranking | 88.2 | 66.92M | 0.1 min |
| ANNOVAR (SNV) | Small-variant annotation | 98.0 | 8.15GB | 151.7 min |
| RankVar | Small-variant prioritization | 96.4 | 631.7M | 1.7 min |
| RankScore | Small-variant scoring | 96.2 | 10.08M | 1.7 min |
| ANNOVAR (SV) | Structural-variant annotation | 91.7 | 507M | 0.4 min |
| SURVIVOR | Structural-variant conversion | 3.07 | 3.074M | 0.0 min |
| PhenoSV | Structural-variant phenotype scoring | 1.150 | 1.15GB | 10.8 min |

### Short-read routes

| Software | Role | CPU usage (reported) | RAM usage (reported) | Time (reported) |
| --- | --- | ---: | ---: | ---: |
| ExpansionHunter | Repeat analysis | 97.7 | 55.33M | 1.3 min |
| Manta | Structural variants | 388.1 | 660.1M | 88.3 min |
| DeepVariant | Default small variants | 1553.5 | 35.90GB | 190.1 min |
| Truvari | Structural consolidation | 70.0 | 106.0M | 0.2 min |
| CNVnator | Copy-number analysis | 99.5 | 13.25GB | 159.6 min |
| HaplotypeCaller | Light small variants | 135.3 | 7.259GB | 978.6 min |

### Long-read routes

| Software | Role | CPU usage (reported) | RAM usage (reported) | Time (reported) |
| --- | --- | ---: | ---: | ---: |
| LongPhase | Combined-analysis phasing | 179.6 | 11.59GB | 24.9 min |
| NanoCaller | Light small variants | 790.1 | 6.509G | 270.3 min |
| Clair3 | Default small variants | 599.3 | 12.59GB | 105.2 min |
| Sniffles2 | Structural variants | 189.8 | 451.2M | 3.9 min |

### Phenotype extraction

| Software | Role | CPU usage (reported) | RAM usage (reported) | Time (reported) |
| --- | --- | ---: | ---: | ---: |
| PhenoTagger | Default extractor | 95.3 | 2.161GB | 6.7 min |
| PhenoGPT2 | Optional extractor | 108.4 | 10.53GB | 1.3 min |

For capacity planning, start with the configured allocations rather than these
observations, then adjust concurrency using measurements from your own data and
infrastructure. Preserve the Nextflow report and trace so later comparisons use
the same measurement definitions.

Return to [Installation and setup](INSTALLATION.md), [Running PipeVar](USAGE.md),
or the [documentation map](../README.md#documentation).
