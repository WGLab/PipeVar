## PipeVar_mito Usage Notes

- Enable MELT with `--melt yes` only for short-read BAM/CRAM runs.
- `--melt yes` is valid with `--mode sv` or when `--mode` is omitted.
- `--melt yes` is rejected for VCF-only input, `--type ont|pacbio`, and `--mode snp`.
- The MELT container must already include `MELT.jar` and bundled hg38/grch38 resources under `/opt/melt/resources`.
- MELT outputs a plain `*_melt.vcf`, which is merged into the short-read SV callset before ANNOVAR, SURVIVOR, PhenoSV, and final prioritization.
- Enable mitochondrial calling with `--mito yes` only for short-read BAM/CRAM runs.
- Mito mode requires a fully indexed reference bundle next to `--ref_fa`: `.fai`, `.dict`, `.amb`, `.ann`, `.bwt`, `.pac`, and `.sa`.
- The mitochondrial prep path keeps the current CRAM decode flow and assumes DRAGEN CRAM input was created against the exact supplied reference bundle.
- The DRAGEN compatibility path uses `gatk RevertSam --RESTORE_HARDCLIPS false` before `bwa mem` realignment.
