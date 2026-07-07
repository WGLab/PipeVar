import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path("/home/beoungle/docker_work/truvari/prepare_shortread_sv_merge.py")
PIPEVAR_MITO = Path("/home/beoungle/PipeVar_mito")


class PrepareShortreadSvMergeTest(unittest.TestCase):
    def run_helper(self, root, vcfs, env=None):
        work_dir = root / "truvari_inputs"
        subprocess.run(
            [
                sys.executable,
                str(SCRIPT),
                "--out-prefix",
                "sample",
                "--reference",
                "/ref.fa",
                "--work-dir",
                str(work_dir),
                "--prepare-only",
                *map(str, vcfs),
            ],
            check=True,
            cwd=root,
            env=env,
        )
        return work_dir

    def test_manta_sampleft_format_is_defined_when_missing_from_header(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            manta = root / "sample_manta.vcf"
            manta.write_text(
                "\n".join(
                    [
                        "##fileformat=VCFv4.2",
                        "##contig=<ID=chr1,length=248956422>",
                        '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Structural variant type">',
                        '##INFO=<ID=END,Number=1,Type=Integer,Description="End coordinate">',
                        '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Structural variant length">',
                        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
                        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample",
                        "chr1\t1427447\tmanta_del_1\tN\t<DEL>\t60\tPASS\tSVTYPE=DEL;END=1428447;SVLEN=-1000\tGT:SampleFT\t0/1:PASS",
                        "",
                    ]
                )
            )

            work_dir = self.run_helper(root, [manta])
            normalized = (work_dir / "manta.normalized.vcf").read_text()

            self.assertIn("##FORMAT=<ID=SampleFT,Number=1,Type=String", normalized)
            self.assertIn("\tGT:SampleFT\t0/1:PASS", normalized)

    def test_prefers_bcftools_view_header_when_available(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            manta = root / "sample_manta.vcf"
            manta.write_text(
                "\n".join(
                    [
                        "##fileformat=VCFv4.2",
                        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample",
                        "chr1\t1427447\tmanta_del_1\tN\t<DEL>\t60\tPASS\tSVTYPE=DEL;END=1428447;SVLEN=-1000\tGT:SampleFT\t0/1:PASS",
                        "",
                    ]
                )
            )
            fake_bin = root / "bin"
            fake_bin.mkdir()
            fake_bcftools = fake_bin / "bcftools"
            fake_bcftools.write_text(
                "\n".join(
                    [
                        "#!/bin/sh",
                        "if [ \"$1\" != \"view\" ] || [ \"$2\" != \"-h\" ]; then exit 2; fi",
                        "printf '%s\\n' '##fileformat=VCFv4.2'",
                        "printf '%s\\n' '##source=fake-bcftools-header'",
                        "printf '%s\\n' '##contig=<ID=chr1,length=248956422>'",
                        "printf '%s\\n' '##FILTER=<ID=LowQual,Description=\"Low quality\">'",
                        "printf '%s\\n' '##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Structural variant type\">'",
                        "printf '%s\\n' '##INFO=<ID=END,Number=1,Type=Integer,Description=\"End coordinate\">'",
                        "printf '%s\\n' '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Structural variant length\">'",
                        "printf '%s\\n' '##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">'",
                        "printf '%s\\n' '##FORMAT=<ID=SampleFT,Number=1,Type=String,Description=\"Sample filter\">'",
                        "printf '%s\\n' '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample'",
                        "",
                    ]
                )
            )
            fake_bcftools.chmod(0o755)
            env = os.environ.copy()
            env["PATH"] = f"{fake_bin}:{env.get('PATH', '')}"

            work_dir = self.run_helper(root, [manta], env=env)
            normalized = (work_dir / "manta.normalized.vcf").read_text()

            self.assertIn("##source=fake-bcftools-header", normalized)
            self.assertIn("##FILTER=<ID=LowQual", normalized)
            self.assertIn("##FORMAT=<ID=SampleFT,Number=1,Type=String", normalized)

    def test_helper_normalizes_only_and_does_not_write_merge_outputs(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            manta = root / "sample_manta.vcf"
            manta.write_text(
                "\n".join(
                    [
                        "##fileformat=VCFv4.2",
                        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
                        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample",
                        "chr1\t1427447\tmanta_del_1\tN\t<DEL>\t60\tPASS\tSVTYPE=DEL;END=1428447;SVLEN=-1000\tGT\t0/1",
                        "",
                    ]
                )
            )

            work_dir = self.run_helper(root, [manta])

            self.assertTrue((work_dir / "manta.normalized.vcf").exists())
            self.assertFalse((root / "sample.shortread_sv.merged.vcf").exists())
            self.assertFalse((root / "sample.shortread_sv.truvari_collapsed.vcf").exists())

    def test_modules_run_merge_and_collapse_directly(self):
        module_paths = [
            PIPEVAR_MITO / "modules/truvari_shortread_sv_merge/main.nf",
            PIPEVAR_MITO / "modules/multi_truvari_shortread_sv_merge/main.nf",
        ]

        helper_source = SCRIPT.read_text()
        self.assertNotIn("truvari collapse", helper_source)
        self.assertNotIn("bcftools merge", helper_source)

        for module_path in module_paths:
            module_source = module_path.read_text()
            self.assertIn("prepare_shortread_sv_merge.py", module_source)
            self.assertIn("bcftools sort", module_source)
            self.assertIn("bcftools merge -m id", module_source)
            self.assertIn("truvari collapse", module_source)


if __name__ == "__main__":
    unittest.main()
