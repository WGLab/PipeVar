## PipeVar_annotated_snv Usage Notes

- Enable pre-annotated SNP input with `--annotated_snv yes`.
- Annotated-SNV mode requires:
  - `--annovar_txt <sample.hg38_multianno.txt>`
  - `--vcf <sample.hg38_multianno.vcf>`
  - one phenotype source: `--note <FILE>` or `--hpo <FILE>`
- Annotated-SNV mode skips runtime ANNOVAR re-annotation.
- If annotated-SNV is combined with BAM/CRAM input and `--type short`, PipeVar_annotated_snv also runs short-read SV, CNV, STR, and optional mito analysis.
- Annotated-SNV mode rejects `--mode sv` and `--target yes`.
- Unified batch manifests now use:
  - `sample,input_kind,phenotype_path,phenotype_format`
  - conditional columns `snv_txt_path,snv_vcf_path,vcf_path,alignment_path,alignment_index_path`
- For `input_kind=annotated_snv`, `alignment_path` and `alignment_index_path` are optional; if present for all rows, the batch run uses the hybrid annotated-SNV + all-NGS path.
- Supported unified `input_kind` values are `annotated_snv`, `vcf_snv`, `vcf_sv`, `bam_ngs`, and `cram_ngs`.

- Enable SCRAMBLE with `--scramble yes` only for short-read BAM/CRAM runs.
- `--scramble yes` is valid with `--mode sv` or when `--mode` is omitted.
- `--scramble yes` is rejected for VCF-only input, `--type ont|pacbio`, and `--mode snp`.
- PipeVar runs SCRAMBLE as one internal short-read step that performs both cluster identification and cluster analysis.
- The shared SCRAMBLE image must already include `cluster_identifier`, `SCRAMble.R`, and `MEI_consensus_seqs.fa`.
- SCRAMBLE outputs a plain `*_scramble.vcf`, which is merged into the short-read SV callset before ANNOVAR, SURVIVOR, PhenoSV, and final prioritization.
- Enable mitochondrial calling with `--mito yes` only for short-read BAM/CRAM runs.
- Mito mode requires a fully indexed reference bundle next to `--ref_fa`: `.fai`, `.dict`, `.amb`, `.ann`, `.bwt`, `.pac`, and `.sa`.
- The mitochondrial prep path keeps the current CRAM decode flow and assumes DRAGEN CRAM input was created against the exact supplied reference bundle.
- The DRAGEN compatibility path uses `gatk RevertSam --RESTORE_HARDCLIPS false` before `bwa mem` realignment.
