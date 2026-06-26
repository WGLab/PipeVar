import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path("/home/beoungle/docker_work/xtea/assets/xtea_candidates_to_vcf.py")


def candidate_line(chrom, left, right, ins_len):
    fields = [
        chrom,
        str(left),
        str(left),
        str(right),
        "2",
        "6",
        "5",
        "4",
        "3",
        "1",
        "0",
    ]
    while len(fields) < 47:
        fields.append("0")
    fields.append(str(ins_len))
    return "\t".join(fields) + "\n"


class XteaCandidatesToVcfTest(unittest.TestCase):
    def run_converter(self, root):
        out_vcf = root / "sample_xtea.vcf"
        log_file = root / "xtea.run.log"
        subprocess.run(
            [
                sys.executable,
                str(SCRIPT),
                "--sample-id",
                "sample",
                "--output",
                str(out_vcf),
                "--reference",
                "/ref.fa",
                "--log",
                str(log_file),
                str(root),
            ],
            check=True,
        )
        return out_vcf.read_text(), log_file.read_text()

    def data_lines(self, text):
        return [line for line in text.splitlines() if line and not line.startswith("#")]

    def test_converts_priority_candidate_files_and_ignores_internal_snp_vcf(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            alu = root / "sample_xtea_work" / "ALU"
            sva = root / "sample_xtea_work" / "SVA"
            alu.mkdir(parents=True)
            sva.mkdir(parents=True)
            (alu / "internal_snp.vcf.gz").write_text("not an MEI VCF\n")
            (alu / "candidate_disc_filtered_cns.txt.high_confident.post_filtering.txt").write_text(
                candidate_line("chr1", 120, 150, 300)
            )
            (sva / "candidate_disc_filtered_cns_post_filtering.txt").write_text(
                candidate_line("chr2", 220, 200, 90)
            )

            vcf, log = self.run_converter(root)
            records = self.data_lines(vcf)

            self.assertEqual(2, len(records))
            self.assertIn("chr1\t120\txtea_1\tN\t<INS:ME:ALU>", records[0])
            self.assertIn("SVTYPE=INS;END=150;SVLEN=300;MEINFO=ALU", records[0])
            self.assertIn("chr2\t200\txtea_2\tN\t<INS:ME:SVA>", records[1])
            self.assertIn("Ignoring xTEA intermediate SNP VCF", log)
            self.assertIn("xTEA gVCF step bypassed; using -f 4883", log)

    def test_empty_candidate_file_emits_header_only_vcf(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            l1 = root / "sample_xtea_work" / "L1"
            l1.mkdir(parents=True)
            (l1 / "candidate_disc_filtered_cns.txt.high_confident.post_filtering.txt").write_text("")

            vcf, log = self.run_converter(root)

            self.assertIn("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample", vcf)
            self.assertEqual([], self.data_lines(vcf))
            self.assertIn("No non-empty xTEA MEI candidate records found", log)


if __name__ == "__main__":
    unittest.main()
