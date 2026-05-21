## PipeVar_mito Usage Notes

- Enable xTEA with `--xtea yes` only for short-read BAM/CRAM runs.
- `--xtea yes` is valid with `--mode sv` or when `--mode` is omitted.
- `--xtea yes` is rejected for VCF-only input, `--type ont|pacbio`, and `--mode snp`.
- PipeVar runs xTEA as one internal short-read step that generates xTEA inputs, runs the local wrapper, and normalizes the per-sample VCF.
- The shared xTEA image must already include the `xtea` command, xTEA scripts, repeat library, and GENCODE GFF3.
- xTEA outputs a plain `*_xtea.vcf`, which is Truvari-merged and deduplicated with Manta and optional CNVnator before ANNOVAR, SURVIVOR, PhenoSV, and final prioritization.
- Enable mitochondrial calling with `--mito yes` only for BAM/CRAM runs.
- `--mito yes` is rejected for VCF-only input, `--mode sv`, and long-read `--light yes` runs.
- Short-read mito uses the existing Mutect2 branch and requires a fully indexed reference bundle next to `--ref_fa`: `.fai`, `.dict`, `.amb`, `.ann`, `.bwt`, `.pac`, and `.sa`.
- Long-read mito uses a mito-specific Clair3 branch and only needs the standard FASTA `.fai` plus BAM/CRAM index.
- The short-read mitochondrial prep path keeps the current CRAM decode flow and assumes DRAGEN CRAM input was created against the exact supplied reference bundle.
- The DRAGEN compatibility path uses `gatk RevertSam --RESTORE_HARDCLIPS false` before `bwa mem` realignment.
- Enable CNVpytor with `--cnvpytor yes` only for long-read BAM/CRAM runs.
- `--cnvpytor yes` is rejected for VCF-only input, `--type short`, and `--mode snp`.
- In `--mode sv`, CNVpytor runs in read-depth-only mode and is merged with Sniffles before `ANNOVAR_SV`.
- In full long-read mode, `--cnvpytor_baf yes` allows CNVpytor to use the branch's existing long-read SNP VCF for SNP/BAF support.
- CNVpytor is experimental for long reads and intended for large CNVs; calls below 100 kb are noisy by default.
- CNVpytor is for whole-genome nuclear long-read analysis only; mitochondrial CNV interpretation is out of scope.
