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
- Long-read BAM/CRAM runs use Sniffles as their structural-variant caller.
- Full ONT and PacBio workflows send the complete Sniffles VCF with `RNAMES` to LongPhase; targeted, de novo, common-SV, and phenotype filters affect the reporting branch, not the phasing context.
- PipeVar retains only the whole-event PhenoSV `Elements=SV` score, imports genes from the curated ANNOVAR `Gene.refGene` field, and intersects both with the phased Sniffles output by exact VCF ID before final prioritization.
- Confirmed trans compound heterozygotes require matching nonmissing `PS`; missing or different phase sets follow `--allow_unphased_comphet`.
