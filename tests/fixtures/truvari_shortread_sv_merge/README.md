# Truvari short-read SV merge fixtures

Tiny VCFs for smoke-testing Manta/CNVnator/xTEA preparation before `bcftools merge` and `truvari collapse --intra`.

- `sample_manta.vcf`: sample VCF with an overlapping DEL and a Manta-only DUP.
- `sample_cnvnator.vcf`: sample-less VCF with a missing ID and overlapping DEL.
- `sample_xtea.vcf`: sample-less MEI/INS-style VCF that should remain distinct from the DEL/DUP calls.

