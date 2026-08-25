# PipeVar parameter reference

This page lists the documented operator-facing parameters for PipeVar 0.5.0.
Defaults are effective values from `main.nf` and `nextflow.config`.

String toggles require explicit values, for example `--mito yes` and
`--GPU yes`; they are not bare flags.

## Inputs, outputs, and route selection

| Parameter | Default | Meaning |
| --- | --- | --- |
| `--bam` | `null` | Single BAM/CRAM input; use `true` to select aligned reads in a basic sample sheet |
| `--vcf` | `null` | Single VCF input; use `true` to select VCF input in a basic sample sheet |
| `--input_csv` | `null` | Basic or typed sample sheet |
| `--annotated_snv` | `no` | Enable prepared ANNOVAR small-variant input |
| `--annovar_txt` | `null` | ANNOVAR hg38 multianno TXT paired with `--vcf` |
| `--ref_fa` | `null` | Reference FASTA |
| `--out_prefix` | `PipeVar` | Single-sample output prefix |
| `--output_directory` | Launch directory | Published-output directory |
| `--type` | `ont` | `short`, `ont`, or `pacbio` for aligned reads |
| `--mode` | `null` | `snp`, `sv`, or omitted for supported combined analysis |
| `--light` | `no` effectively | Select documented light behavior with `yes` |
| `--help` | `false` | Print compact command help and exit |

See the [input guide](INPUTS.md) for required file combinations.

## Phenotype and prioritization

| Parameter | Default | Meaning |
| --- | --- | --- |
| `--note` | `null` | Clinical note; `no` means basic-sheet `note_path` contains HPO terms |
| `--hpo` | `null` | HPO-term file |
| `--inheritance_mode` | `ml` | `ml`, `omim`, or `gnomad`/LOEUF fallback behavior |
| `--gnomad_af_ad` | `0.001` | Dominant small-variant frequency ceiling |
| `--gnomad_af_ar` | `0.01` | Recessive and upstream small-variant frequency ceiling |
| `--include_clinvar_report` | `yes` | Include candidates supported only by accepted ClinVar evidence |
| `--allow_unphased_comphet` | `no` | Allow unresolved compound-heterozygous pairs |
| `--prioritize_sv_only` | `no` | Limit combined final prioritization to structural-variant evidence |
| `--rankscore` | `0.50` | Minimum RankScore value |
| `--rankscore_softwares` | All built-ins | Comma-separated RankScore software subset |
| `--rankvar` | `0.05` | Minimum RankVar score |
| `--phenosv_score` | `0.50` | Minimum PhenoSV score |
| `--gq` | `20` | Minimum genotype quality |
| `--ad` | `15` | Minimum allele depth |
| `--nanocaller_dp` | `20` | NanoCaller depth fallback when GQ is absent |
| `--phen2gene_filter` | `500` | Number of top Phen2Gene genes used for targeted regions |
| `--target` | `no` effectively | Enable phenotype-derived targeted analysis with `yes` |
| `--gene` | `null` | Comma-separated genes or one-gene-per-line file |

The dominant frequency ceiling must not exceed the recessive ceiling.

## Callers and repeat analysis

| Parameter | Default | Meaning |
| --- | --- | --- |
| `--cnvnator` | `yes` | Add CNVnator in supported short-read structural/combined routes |
| `--cnvnator_bin_size` | `100` | CNVnator read-depth bin size |
| `--genome` | `hg38` | Bundled ExpansionHunter catalog naming (`hg38` or `grch38`) |
| `--expansionhunter_variant_catalog` | `null` | Override the bundled ExpansionHunter catalog |
| `--truvari_shortread_refdist` | `1000` | Short-read structural-variant merge reference distance |
| `--truvari_shortread_pctseq` | `0` | Merge sequence-similarity threshold |
| `--truvari_shortread_pctsize` | `0.5` | Merge size-similarity threshold |
| `--truvari_shortread_pctovl` | `0.5` | Merge overlap threshold |
| `--truvari_shortread_sizemin` | `50` | Minimum variant size for collapsing |
| `--truvari_shortread_keep` | `first` | Record-selection rule during collapsing |

## Common structural-variant filtering

| Parameter | Default | Meaning |
| --- | --- | --- |
| `--common_sv_filter` | `yes` | Enable common-variant annotation and inheritance-aware filtering |
| `--common_sv_af_ad` | `0.005` effective | Dominant and unknown-model frequency ceiling |
| `--common_sv_af_ar` | `0.01` effective | Recessive and upstream frequency ceiling |
| `--common_sv_af` | `null` | Deprecated alias that sets both ceilings when used alone |
| `--common_sv_reciprocal_overlap` | `0.5` | Reciprocal overlap for interval matches |
| `--common_sv_distance` | `1000` | Breakpoint fallback distance |
| `--common_sv_ins_distance` | `500` | Insertion position window |
| `--common_sv_ins_identity` | `0.5` | Insertion sequence-identity threshold |

`--common_sv_af` cannot be combined with either split frequency option.

## Mitochondrial analysis

| Parameter | Default | Meaning |
| --- | --- | --- |
| `--mito` | `no` | Enable mitochondrial analysis with `yes` |
| `--mito_contig` | `chrM` | Preferred mitochondrial contig alias |
| `--mito_min_vaf` | `0.01` | Prioritization variant-allele-fraction floor |
| `--mito_min_depth` | `50` | Prioritization depth floor |
| `--mito_min_alt_reads` | `5` | Prioritization alternate-read floor |
| `--mito_gui_min_af` | `0.5` | HTML-report allele-fraction threshold |
| `--mito_gui_min_apogee2` | `0.5` | HTML-report APOGEE2 threshold |
| `--mito_gui_min_mitotip` | `12.66` | HTML-report MitoTip threshold |
| `--hmtvar_data` | `null` | Optional HmtVar data state/input |

## Runtime and containers

| Parameter | Default | Meaning |
| --- | --- | --- |
| `--annovar_host_path` | `${projectDir}/annovar` | ANNOVAR host directory mounted at `/annovar` |
| `--phenosv_host_path` | `${projectDir}/PhenoSV_model` | PhenoSV host directory mounted at `/PhenoSV/train_data` |
| `--GPU` | `no` | Enable GPU options; preserve the uppercase parameter name |
| `--gpu_backend` | `singularity` | `singularity` or `docker`; profiles set this automatically |
| `--gpu_cpus` | `16` | CPU allocation for GPU-backed tasks |
| `--gpu_cluster_options` | `--gres=gpu:1` | Scheduler options for GPU tasks |
| `--deepvariant_max_forks` | Effectively unlimited | Optional DeepVariant concurrency limit |

## Phenotype extraction

| Parameter | Default | Meaning |
| --- | --- | --- |
| `--phenotype_extractor` | `phenotagger` | `phenotagger` or `phenogpt2` |
| `--phenogpt2_batch_size` | `1` | Supported inference batch size |
| `--phenogpt2_chunk_batch_size` | `1` | Supported chunk batch size |
| `--phenogpt2_wc` | `0` | Word-count chunking; supported value disables it |
| `--phenogpt2_attn_implementation` | `eager` | Attention implementation |
| `--phenogpt2_negation` | `no` | Enable negation and embedding verification |
| `--phenogpt2_max_forks` | `1` | Maximum concurrent PhenoGPT2 tasks |
| `--phenogpt2_model_host_path` | `null` | Required external base-model directory |
| `--phenogpt2_negation_model_host_path` | `null` | Required with negation processing |
| `--phenogpt2_embedding_model_host_path` | `null` | Required with negation processing |
| `--phenogpt2_cache_host_path` | `null` | Optional pre-created writable cache |

See [Installation and setup](INSTALLATION.md#phenogpt2-and-gpu-execution) for
mount and cluster requirements.

Return to [Running PipeVar](USAGE.md) or the
[documentation map](../README.md#documentation).
