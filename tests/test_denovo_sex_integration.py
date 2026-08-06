import csv
import importlib.util
import os
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]


def load_module(name, path):
    module_dir = str(path.parent)
    if module_dir not in sys.path:
        sys.path.insert(0, module_dir)
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


DENOVO = load_module(
    "pipevar_denovo_filter",
    ROOT / "docker_work/denovo_snv_sv_filter/filter_denovo_snv_sv.py",
)
PRIO = load_module(
    "pipevar_sex_prio",
    ROOT / "docker_work/longphase/prio_gene_only.py",
)


class CombinedDenovoSexTests(unittest.TestCase):
    def test_inherited_sex_linked_calls_are_removed_before_male_promotion(self):
        with tempfile.TemporaryDirectory() as tmp:
            work = Path(tmp)
            (work / "trio.csv").write_text(
                "sample,family_id,role,vcf_sample,sex\n"
                "child,F1,proband,KID,male\n"
                "mom,F1,mother,MOM,female\n"
                "dad,F1,father,DAD,male\n"
            )
            header = (
                "##fileformat=VCFv4.2\n"
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tKID\tMOM\tDAD\n"
            )
            inherited_x = "Gene=GENEX;SOURCE=rankvar;MODEL=XLR;RANK_VAR=0.80"
            novel_x = "Gene=GENEX;SOURCE=rankvar;MODEL=XLR;RANK_VAR=0.91"
            inherited_y = "Gene=GENEY;SOURCE=rankvar;MODEL=Recessive;RANK_VAR=0.75"
            novel_y = "Gene=GENEY;SOURCE=rankvar;MODEL=Recessive;RANK_VAR=0.88"
            (work / "child.vcf").write_text(
                header
                + f"chrX\t100\t.\tA\tG\t.\tPASS\t{inherited_x}\tGT\t1\t0\t0\n"
                + f"chrX\t200\t.\tC\tT\t.\tPASS\t{novel_x}\tGT\t1\t0\t0\n"
                + f"chrY\t300\t.\tG\tA\t.\tPASS\t{inherited_y}\tGT\t1\t0\t0\n"
                + f"chrY\t400\t.\tT\tC\t.\tPASS\t{novel_y}\tGT\t1\t0\t0\n"
            )
            (work / "mom.vcf").write_text(
                header + f"X\t100\t.\tA\tG\t.\tPASS\t{inherited_x}\tGT\t0\t1\t0\n"
            )
            (work / "dad.vcf").write_text(
                header + f"Y\t300\t.\tG\tA\t.\tPASS\t{inherited_y}\tGT\t0\t0\t1\n"
            )

            previous = Path.cwd()
            try:
                os.chdir(work)
                DENOVO.main(
                    [
                        "--mode", "snv",
                        "--pedigree-csv", "trio.csv",
                        "--snv-vcf-input", "child", "child.vcf",
                        "--snv-vcf-input", "mom", "mom.vcf",
                        "--snv-vcf-input", "dad", "dad.vcf",
                    ]
                )
            finally:
                os.chdir(previous)

            output = work / "child.denovo.snv.vcf"
            text = output.read_text()
            self.assertNotIn("\t100\t", text)
            self.assertNotIn("\t300\t", text)
            self.assertIn("\t200\t", text)
            self.assertIn("\t400\t", text)

            _, variants = PRIO.parse_vcf(output)
            categories = {
                scenario["cat"]
                for gene_variants in (
                    [variant for variant in variants if gene in variant.genes]
                    for gene in ("GENEX", "GENEY")
                )
                for scenario in PRIO.get_all_gene_scenarios(gene_variants, sex="male")
            }
            self.assertIn("RankVar_AR_Hemizygous", categories)
            with (work / "denovo.snv.summary.tsv").open() as handle:
                summary = next(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(summary["variants_removed"], "2")
            self.assertEqual(summary["variants_kept"], "2")


if __name__ == "__main__":
    unittest.main()
