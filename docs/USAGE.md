## PipeVar_mito Usage Notes

- Enable SCRAMBLE with `--scramble yes` only for short-read BAM/CRAM runs.
- `--scramble yes` is valid with `--mode sv` or when `--mode` is omitted.
- `--scramble yes` is rejected for VCF-only input, `--type ont|pacbio`, and `--mode snp`.
- PipeVar runs SCRAMBLE as direct `clusteridentifier` and `clusteranalysis` stages inside the existing short-read subworkflows.
- The shared SCRAMBLE image must already include `cluster_identifier`, `SCRAMble.R`, and `MEI_consensus_seqs.fa`.
- SCRAMBLE outputs a plain `*_scramble.vcf`, which is merged into the short-read SV callset before ANNOVAR, SURVIVOR, PhenoSV, and final prioritization.
- Enable mitochondrial calling with `--mito yes` only for short-read BAM/CRAM runs.
- Mito mode requires a fully indexed reference bundle next to `--ref_fa`: `.fai`, `.dict`, `.amb`, `.ann`, `.bwt`, `.pac`, and `.sa`.
- The mitochondrial prep path keeps the current CRAM decode flow and assumes DRAGEN CRAM input was created against the exact supplied reference bundle.
- The DRAGEN compatibility path uses `gatk RevertSam --RESTORE_HARDCLIPS false` before `bwa mem` realignment.
