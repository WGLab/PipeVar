import csv
import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path


MODULE_PATH = Path(__file__).parents[1] / "scripts" / "compare_gene_inheritance.py"
SPEC = importlib.util.spec_from_file_location("compare_gene_inheritance", MODULE_PATH)
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


class InheritanceComparisonTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.pipevar = self.root / "pipevar"
        self.exomiser = self.root / "exomiser"
        self.pipevar.mkdir()
        self.exomiser.mkdir()

    def tearDown(self):
        self.temp.cleanup()

    def write(self, path, text):
        path.write_text(text, encoding="utf-8")
        return path

    def test_truth_normalizes_models_and_merges_distinct_pairs(self):
        truth = self.write(
            self.root / "truth.tsv",
            "sample_id\tgenes\tinheritance\n"
            "S1\tGENE1\tautosomal recessive\n"
            "S1\tGENE1\tdominant\n"
            "S1\tGENE2;GENE3\tx-linked\n"
            "S2\tGENE4\tnot specified\n",
        )
        samples = MODULE.read_truth_tsv(truth)
        self.assertEqual(
            (
                MODULE.TruthCandidate("GENE1", "AR"),
                MODULE.TruthCandidate("GENE1", "AD"),
                MODULE.TruthCandidate("GENE2", "XL"),
                MODULE.TruthCandidate("GENE3", "XL"),
            ),
            samples[0].candidates,
        )
        self.assertFalse(MODULE.SampleResult(
            samples[1], Path(), Path(), MODULE.ParsedResult("ok", ()),
            MODULE.ParsedResult("ok", ()), MODULE.ToolMatch(), MODULE.ToolMatch()
        ).inheritance_evaluable)

    def test_truth_rejects_duplicate_gene_model_pair(self):
        truth = self.write(
            self.root / "duplicate.tsv",
            "sample_id\tgenes\tinheritance\nS1\tGENE1\tAR\nS1\tgene1\trecessive\n",
        )
        with self.assertRaisesRegex(MODULE.ComparisonError, "duplicate GENE1/AR"):
            MODULE.read_truth_tsv(truth)

    def test_pipevar_prefers_explicit_xlinked_model_from_vcf_tail(self):
        path = self.write(
            self.root / "pipevar.tsv",
            "#RANK\tGENE\tPRIORITY_TIER\tVCF_LINE\n"
            "1\tGENEX\tRankVar_AR_Hemizygous\t"
            "chrX\t100\t.\tA\tG\t.\tPASS\tGene=GENEX;MODEL=XLR\tGT\t0/1\n"
            "2\tGENEA\tRankVar_AD\t"
            "chr1\t200\t.\tC\tT\t.\tPASS\tGene=GENEA;MODEL=Dominant\tGT\t0/1\n",
        )
        rankings = MODULE.read_pipevar(path)
        self.assertEqual([("GENEX", "XLR", 1), ("GENEA", "AD", 2)], [
            (item.gene, item.inheritance, item.rank) for item in rankings
        ])

    def test_pipevar_recovers_multiple_encoded_or_legacy_models(self):
        encoded = MODULE.parse_info_models(
            "chr1\t1\t.\tA\tG\t.\tPASS\tGene=G;MODEL=Dominant%3bRecessive"
        )
        legacy = MODULE.parse_info_models(
            "chr1\t1\t.\tA\tG\t.\tPASS\tGene=G;MODEL=Dominant;Recessive"
        )
        self.assertEqual(("AD", "AR"), encoded)
        self.assertEqual(("AD", "AR"), legacy)

    def test_exomiser_accepts_duplicate_ranks_and_x_recessive(self):
        path = self.write(
            self.root / "genes.tsv",
            "#RANK\tGENE_SYMBOL\tMOI\n"
            "1\tGENE1\tAD\n1\tGENE1\tAR\n2\tGENEX\tX_RECESSIVE\n",
        )
        rankings = MODULE.read_exomiser(path)
        self.assertEqual(["AD", "AR", "XLR"], [item.inheritance for item in rankings])

    def test_generic_xlinked_truth_matches_xld_or_xlr(self):
        self.assertTrue(MODULE.inheritance_matches("XL", "XLD"))
        self.assertTrue(MODULE.inheritance_matches("XL", "XLR"))
        self.assertFalse(MODULE.inheritance_matches("XL", "AD"))

    def test_end_to_end_reports_exact_and_gene_only_accuracy(self):
        truth = self.write(
            self.root / "truth.tsv",
            "sample_id\tgenes\tinheritance\nS1\tGENE1\tAR\nS2\tGENEX\tXL\n",
        )
        self.write(
            self.pipevar / "S1.norm.alt.prio_gene.vcf",
            "#RANK\tGENE\tPRIORITY_TIER\tVCF_LINE\n"
            "1\tGENE1\tRankVar_AD\tchr1\t1\t.\tA\tG\t.\tPASS\tGene=GENE1;MODEL=Dominant\tGT\t0/1\n",
        )
        self.write(
            self.exomiser / "S1.exomiser.input.vcf.gz-results.genes.tsv",
            "#RANK\tGENE_SYMBOL\tMOI\n1\tGENE1\tAR\n",
        )
        self.write(
            self.pipevar / "S2.norm.alt.prio_gene.vcf",
            "#RANK\tGENE\tPRIORITY_TIER\tVCF_LINE\n"
            "1\tGENEX\tRankVar_AR_Hemizygous\tchrX\t2\t.\tC\tT\t.\tPASS\tGene=GENEX;MODEL=XLR\tGT\t0/1\n",
        )
        self.write(
            self.exomiser / "S2.exomiser.input.vcf.gz-results.genes.tsv",
            "#RANK\tGENE_SYMBOL\tMOI\n1\tGENEX\tXR\n",
        )
        exit_code = MODULE.main([
            "--truth-tsv", str(truth), "--pipevar-dir", str(self.pipevar),
            "--exomiser-dir", str(self.exomiser),
            "--output-prefix", str(self.root / "comparison"),
        ])
        self.assertEqual(0, exit_code)
        with (self.root / "comparison.inheritance.summary.tsv").open(newline="") as handle:
            summary = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(["1", "5", "10", "20"], [row["cutoff"] for row in summary])
        self.assertEqual("2", summary[0]["paired_inheritance_evaluable_samples"])
        self.assertEqual("1", summary[0]["pipevar_correct_gene_and_inheritance_hits"])
        self.assertEqual("2", summary[0]["pipevar_gene_only_hits"])
        self.assertEqual("2", summary[0]["exomiser_correct_gene_and_inheritance_hits"])
        with (self.root / "comparison.inheritance.details.tsv").open(newline="") as handle:
            details = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual("AD", details[0]["pipevar_gene_predicted_models"])
        self.assertEqual("", details[0]["pipevar_exact_rank"])
        self.assertEqual("1", details[0]["pipevar_gene_only_rank"])


if __name__ == "__main__":
    unittest.main()
