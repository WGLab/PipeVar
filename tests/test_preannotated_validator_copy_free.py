import gzip
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
VALIDATOR = (
    ROOT.parent
    / "docker_work/rankscore/validate_preannotated_annovar_pair.py"
)

TXT_HEADER = (
    "Chr\tStart\tEnd\tRef\tAlt\tFunc.refGene\tGene.refGene\t"
    "Otherinfo4\tOtherinfo5\tOtherinfo7\tOtherinfo8\tOtherinfo12\tOtherinfo13\n"
)


def txt_record(pos: int = 100, format_value: str = "GT:DP", sample_value: str = "0/1:30") -> str:
    return (
        f"1\t{pos}\t{pos}\tA\tG\texonic\tGENE1\t"
        f"1\t{pos}\tA\tG\t{format_value}\t{sample_value}\n"
    )


def vcf_content(
    pos: int = 100,
    format_value: str = "GT:DP",
    sample_value: str = "0/1:30",
    sample_header: str = "SAMPLE",
    annovar_header: bool = True,
) -> str:
    header = "##fileformat=VCFv4.2\n"
    if annovar_header:
        header += (
            '##INFO=<ID=Func.refGene,Number=.,Type=String,'
            'Description="ANNOVAR functional category">\n'
        )
    return (
        header
        + f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample_header}\n"
        + f"1\t{pos}\t.\tA\tG\t.\tPASS\tFunc.refGene=exonic\t"
        + f"{format_value}\t{sample_value}\n"
    )


class PreannotatedValidatorCopyFreeTests(unittest.TestCase):
    def run_validator(self, txt_path: Path, vcf_path: Path) -> subprocess.CompletedProcess[str]:
        return subprocess.run(
            [
                "python3",
                str(VALIDATOR),
                "--sample",
                "case1",
                "--annovar-txt",
                str(txt_path),
                "--annovar-vcf",
                str(vcf_path),
            ],
            text=True,
            capture_output=True,
            check=False,
        )

    def write_valid_pair(self, directory: Path) -> tuple[Path, Path]:
        txt_path = directory / "case1.hg38_multianno.txt"
        vcf_path = directory / "case1.hg38_multianno.vcf"
        txt_path.write_text(TXT_HEADER + txt_record(), encoding="utf-8")
        vcf_path.write_text(vcf_content(), encoding="utf-8")
        return txt_path, vcf_path

    def test_valid_pair_is_unchanged_and_creates_no_files(self):
        with tempfile.TemporaryDirectory() as tmp_name:
            tmp = Path(tmp_name)
            txt_path, vcf_path = self.write_valid_pair(tmp)
            before_names = {path.name for path in tmp.iterdir()}
            before_txt = txt_path.read_bytes()
            before_vcf = vcf_path.read_bytes()

            result = self.run_validator(txt_path, vcf_path)

            self.assertEqual(0, result.returncode, result.stderr)
            self.assertEqual(before_names, {path.name for path in tmp.iterdir()})
            self.assertEqual(before_txt, txt_path.read_bytes())
            self.assertEqual(before_vcf, vcf_path.read_bytes())

    def test_invalid_pairs_still_fail(self):
        cases = {
            "missing_txt_columns": (
                "Chr\tStart\tRef\tAlt\n1\t100\tA\tG\n",
                vcf_content(),
                ".vcf",
                "missing required ANNOVAR columns",
            ),
            "missing_gt": (
                TXT_HEADER + txt_record(format_value="DP", sample_value="30"),
                vcf_content(),
                ".vcf",
                "missing required GT in FORMAT",
            ),
            "non_annovar_vcf": (
                TXT_HEADER + txt_record(),
                vcf_content(annovar_header=False),
                ".vcf",
                "does not look like an ANNOVAR multianno VCF",
            ),
            "multiple_samples": (
                TXT_HEADER + txt_record(),
                vcf_content(sample_header="SAMPLE1\tSAMPLE2").replace("0/1:30\n", "0/1:30\t0/1:30\n"),
                ".vcf",
                "must contain exactly one sample column",
            ),
            "compressed_vcf": (
                TXT_HEADER + txt_record(),
                vcf_content(),
                ".vcf.gz",
                "requires an uncompressed ANNOVAR multianno VCF",
            ),
            "mismatched_pair": (
                TXT_HEADER + txt_record(pos=100),
                vcf_content(pos=101),
                ".vcf",
                "does not match the supplied ANNOVAR TXT",
            ),
        }

        for name, (txt_text, vcf_text, vcf_suffix, expected_error) in cases.items():
            with self.subTest(name=name), tempfile.TemporaryDirectory() as tmp_name:
                tmp = Path(tmp_name)
                txt_path = tmp / "case1.hg38_multianno.txt"
                vcf_path = tmp / f"case1.hg38_multianno{vcf_suffix}"
                txt_path.write_text(txt_text, encoding="utf-8")
                if vcf_suffix == ".vcf.gz":
                    with gzip.open(vcf_path, "wt", encoding="utf-8") as handle:
                        handle.write(vcf_text)
                else:
                    vcf_path.write_text(vcf_text, encoding="utf-8")

                result = self.run_validator(txt_path, vcf_path)

                self.assertNotEqual(0, result.returncode)
                self.assertIn(expected_error, result.stderr)


if __name__ == "__main__":
    unittest.main()
