import subprocess
import tempfile
import unittest
from pathlib import Path


IMAGE = "beoungl/docker_test:common_sv_filter_0.3"


COMMON_VCF = """##fileformat=VCFv4.2
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="SV type">
##INFO=<ID=END,Number=1,Type=Integer,Description="End">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="SV length">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">
##INFO=<ID=SUPPORT,Number=1,Type=Integer,Description="Support">
##INFO=<ID=FUNC,Number=1,Type=String,Description="Function">
##INFO=<ID=GENE,Number=1,Type=String,Description="Gene">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t1000\tCOMMON_DEL\tN\t<DEL>\t.\tMixed\tSVTYPE=DEL;END=2000;SVLEN=-1001;AF=0.02;SUPPORT=12;FUNC=intronic;GENE=GENE1
chr1\t5000\tCOMMON_INS\tA\tATTTT\t.\tPASS\tSVTYPE=INS;END=5000;SVLEN=4;AF=0.05;SUPPORT=20;FUNC=intronic;GENE=GENE1
chr1\t7000\tCOMMON_LOW_AF\tN\t<DUP>\t.\tPASS\tSVTYPE=DUP;END=7800;SVLEN=801;AF=0.005;SUPPORT=2;FUNC=intronic;GENE=GENE1
chr2\t9000\tCOMMON_CHR2\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=9400;SVLEN=-401;AF=0.03;SUPPORT=8;FUNC=intronic;GENE=GENE1
chr1\t12000\tCOMMON_BND\tN\tN[chr2:13000[\t.\tPASS\tSVTYPE=BND;END=12000;SVLEN=1;AF=0.10;SUPPORT=9;FUNC=intronic;GENE=GENE1
"""


def patient_vcf(body: str) -> str:
    return """##fileformat=VCFv4.2
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="SV type">
##INFO=<ID=END,Number=1,Type=Integer,Description="End">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="SV length">
##INFO=<ID=FUNC,Number=1,Type=String,Description="Function">
##INFO=<ID=GENE,Number=1,Type=String,Description="Gene">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE
""" + body


class FilterCommonSVsImageSmokeTest(unittest.TestCase):
    def run_filter(self, sample_text: str):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            common = tmp / "common.vcf"
            sample = tmp / "sample.vcf"
            kept = tmp / "kept.vcf"
            removed = tmp / "removed.vcf"
            summary = tmp / "summary.tsv"
            common.write_text(COMMON_VCF)
            sample.write_text(sample_text)
            subprocess.run(
                [
                    "docker",
                    "run",
                    "--rm",
                    "-v",
                    f"{tmp}:/work",
                    IMAGE,
                    "python3",
                    "/usr/local/bin/filter_common_svs.py",
                    "--input-vcf",
                    "/work/sample.vcf",
                    "--common-vcf",
                    "/work/common.vcf",
                    "--output-vcf",
                    "/work/kept.vcf",
                    "--removed-vcf",
                    "/work/removed.vcf",
                    "--summary-tsv",
                    "/work/summary.tsv",
                    "--common-af",
                    "0.01",
                ],
                check=True,
            )
            return kept.read_text(), removed.read_text(), summary.read_text()

    def data_lines(self, text: str):
        return [line for line in text.splitlines() if line and not line.startswith("#")]

    def test_common_del_removed_and_rare_del_retained(self):
        sample = patient_vcf(
            "chr1\t1100\tPAT_COMMON_DEL\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=1900;SVLEN=-801;FUNC=intronic;GENE=GENE1\tGT\t0/1\n"
            "chr1\t3000\tPAT_RARE_DEL\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=3300;SVLEN=-301;FUNC=intronic;GENE=GENE1\tGT\t0/1\n"
        )
        kept, removed, summary = self.run_filter(sample)
        self.assertIn("PAT_RARE_DEL", kept)
        self.assertNotIn("PAT_COMMON_DEL", kept)
        self.assertIn("PAT_COMMON_DEL", removed)
        self.assertIn("COMMON_SV_ID=COMMON_DEL", removed)
        self.assertIn("removed_DEL\t1", summary)

    def test_reciprocal_overlap_threshold_boundary(self):
        sample = patient_vcf(
            "chr1\t2501\tBELOW_RO\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=3501;SVLEN=-1001;FUNC=intronic;GENE=GENE1\tGT\t0/1\n"
            "chr1\t1400\tABOVE_RO\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=2200;SVLEN=-801;FUNC=intronic;GENE=GENE1\tGT\t0/1\n"
        )
        kept, removed, _summary = self.run_filter(sample)
        self.assertIn("BELOW_RO", kept)
        self.assertIn("ABOVE_RO", removed)

    def test_common_insertion_sequence_removed_but_different_sequence_retained(self):
        sample = patient_vcf(
            "chr1\t5000\tINS_SAME\tA\tATTTT\t.\tPASS\tSVTYPE=INS;END=5000;SVLEN=4;FUNC=intronic;GENE=GENE1\tGT\t0/1\n"
            "chr1\t5000\tINS_DIFF\tA\tAGGGG\t.\tPASS\tSVTYPE=INS;END=5000;SVLEN=4;FUNC=intronic;GENE=GENE1\tGT\t0/1\n"
        )
        kept, removed, _summary = self.run_filter(sample)
        self.assertIn("INS_DIFF", kept)
        self.assertNotIn("INS_SAME", kept)
        self.assertIn("INS_SAME", removed)
        self.assertIn("COMMON_SV_MATCH=INS_SEQUENCE_IDENTITY", removed)

    def test_chromosome_normalization_and_low_af_skip(self):
        sample = patient_vcf(
            "2\t9050\tNO_CHR_DEL\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=9350;SVLEN=-301;FUNC=intronic;GENE=GENE1\tGT\t0/1\n"
            "chr1\t7100\tLOW_AF_DUP\tN\t<DUP>\t.\tPASS\tSVTYPE=DUP;END=7700;SVLEN=601;FUNC=intronic;GENE=GENE1\tGT\t0/1\n"
        )
        kept, removed, _summary = self.run_filter(sample)
        self.assertIn("NO_CHR_DEL", removed)
        self.assertIn("LOW_AF_DUP", kept)

    def test_bnd_retained_and_header_only_kept_vcf_valid(self):
        sample = patient_vcf(
            "chr1\t12000\tPAT_BND\tN\tN[chr2:13000[\t.\tPASS\tSVTYPE=BND;END=12000;SVLEN=1;FUNC=intronic;GENE=GENE1\tGT\t0/1\n"
            "chr1\t1100\tONLY_COMMON\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=1900;SVLEN=-801;FUNC=intronic;GENE=GENE1\tGT\t0/1\n"
        )
        kept, removed, _summary = self.run_filter(sample)
        self.assertIn("#CHROM", kept)
        self.assertIn("PAT_BND", kept)
        self.assertIn("ONLY_COMMON", removed)

    def test_all_removed_kept_vcf_is_header_only(self):
        sample = patient_vcf(
            "chr1\t1100\tONLY_COMMON\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=1900;SVLEN=-801;FUNC=intronic;GENE=GENE1\tGT\t0/1\n"
        )
        kept, removed, summary = self.run_filter(sample)
        self.assertIn("#CHROM", kept)
        self.assertEqual([], self.data_lines(kept))
        self.assertIn("ONLY_COMMON", removed)
        self.assertIn("kept_records\t0", summary)


if __name__ == "__main__":
    unittest.main()
