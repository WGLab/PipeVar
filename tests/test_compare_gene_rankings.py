import argparse
import csv
import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path


MODULE_PATH = Path(__file__).parents[1] / "scripts" / "compare_gene_rankings.py"
SPEC = importlib.util.spec_from_file_location("compare_gene_rankings", MODULE_PATH)
BENCHMARK = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = BENCHMARK
SPEC.loader.exec_module(BENCHMARK)


class CompareGeneRankingsTests(unittest.TestCase):
    def setUp(self):
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary_directory.name)
        self.pipevar_dir = self.root / "pipevar"
        self.exomiser_dir = self.root / "exomiser"
        self.pipevar_dir.mkdir()
        self.exomiser_dir.mkdir()

    def tearDown(self):
        self.temporary_directory.cleanup()

    def write(self, path, text):
        path.write_text(text, encoding="utf-8")
        return path

    def test_truth_supports_comma_and_semicolon_lists(self):
        truth = self.write(
            self.root / "truth.tsv",
            "sample_id\tgenes\nS1\t gene1,GENE2; gene1 \n",
        )
        samples = BENCHMARK.read_truth_tsv(truth)
        self.assertEqual([BENCHMARK.TruthSample("S1", ("GENE1", "GENE2"))], samples)

    def test_truth_merges_repeated_samples_with_distinct_genes(self):
        truth = self.write(
            self.root / "repeated_sample.tsv",
            "sample_id\tgenes\nS1\tGENE1\nS2\tGENE4\nS1\tGENE2;GENE3\n",
        )
        samples = BENCHMARK.read_truth_tsv(truth)
        self.assertEqual(
            [
                BENCHMARK.TruthSample("S1", ("GENE1", "GENE2", "GENE3")),
                BENCHMARK.TruthSample("S2", ("GENE4",)),
            ],
            samples,
        )

    def test_truth_rejects_repeated_gene_for_same_sample(self):
        duplicate = self.write(
            self.root / "duplicate_gene.tsv",
            "sample_id\tgenes\nS1\tGENE1\nS1\tgene1,GENE2\n",
        )
        with self.assertRaisesRegex(
            BENCHMARK.BenchmarkError,
            "duplicate gene 'GENE1' for sample_id 'S1'",
        ):
            BENCHMARK.read_truth_tsv(duplicate)

    def test_truth_rejects_empty_genes(self):
        empty = self.write(
            self.root / "empty_gene.tsv", "sample_id\tgenes\nS1\t ; , \n"
        )
        with self.assertRaisesRegex(BENCHMARK.BenchmarkError, "no gene symbols"):
            BENCHMARK.read_truth_tsv(empty)

    def test_pipevar_parser_handles_current_header_and_embedded_vcf_tabs(self):
        path = self.write(
            self.root / "sample.prio_gene.vcf",
            "#RANK  \tGENE      \tPRIORITY_TIER\tMAX_SCORE\tRANKVAR\t"
            "RANKSCORE\tPHENOSV\tVCF_LINE\n"
            "2\tgene2\ttier\t0.8\t.\t.\t.\tchr2\t200\t.\tC\tT\n"
            "1\tGENE1\ttier\t0.9\t.\t.\t.\tchr1\t100\t.\tA\tG\n",
        )
        rankings = BENCHMARK.read_pipevar_ranking(path)
        self.assertEqual([("GENE1", 1), ("GENE2", 2)], [(r.gene, r.rank) for r in rankings])

    def test_pipevar_parser_handles_legacy_header_and_header_only_result(self):
        legacy = self.write(
            self.root / "legacy.prio_gene.vcf",
            "#RANK\tGENE\tPRIORITY_TIER\tMAX_SCORE\tRANKVAR\tRANKSCORE\tVCF_LINE\n"
            "1\tGENE1\tRankVar_AD\t0.9\t0.9\t.\tchr1\t100\t.\tA\tG\n",
        )
        self.assertEqual("GENE1", BENCHMARK.read_pipevar_ranking(legacy)[0].gene)
        header_only = self.write(
            self.root / "empty.prio_gene.vcf", "#RANK\tGENE\tVCF_LINE\n"
        )
        self.assertEqual((), BENCHMARK.read_pipevar_ranking(header_only))

    def test_pipevar_parser_rejects_duplicate_and_malformed_ranks(self):
        duplicate = self.write(
            self.root / "duplicate_rank.vcf",
            "#RANK\tGENE\n1\tGENE1\n1\tGENE2\n",
        )
        with self.assertRaisesRegex(BENCHMARK.RankingError, "duplicate rank"):
            BENCHMARK.read_pipevar_ranking(duplicate)
        malformed = self.write(
            self.root / "malformed_rank.vcf", "#RANK\tGENE\nfirst\tGENE1\n"
        )
        with self.assertRaisesRegex(BENCHMARK.RankingError, "non-integer rank"):
            BENCHMARK.read_pipevar_ranking(malformed)

    def test_exomiser_deduplicates_moi_rows_and_preserves_raw_rank(self):
        path = self.write(
            self.root / "results.genes.tsv",
            "#RANK\tID\tGENE_SYMBOL\tEXOMISER_GENE_COMBINED_SCORE\n"
            "3\tB_AR\tGENEB\t0.7\n"
            "1\tA_AD\tGENEA\t0.9\n"
            "1\tA_AR\tgenea\t0.8\n"
            "1\tA_XLR\tGENEA\t0.8\n",
        )
        rankings = BENCHMARK.read_exomiser_ranking(path)
        self.assertEqual(
            [("GENEA", 1, 1), ("GENEB", 3, 3)],
            [(r.gene, r.rank, r.raw_rank) for r in rankings],
        )

    def test_exomiser_topk_uses_native_rank_without_gene_rank_compression(self):
        path = self.write(
            self.root / "native_rank.genes.tsv",
            "#RANK\tGENE_SYMBOL\n1\tOTHER\n21\tCAUSAL\n",
        )
        rankings = BENCHMARK.read_exomiser_ranking(path)
        match = BENCHMARK.find_best_match(("CAUSAL",), rankings)
        self.assertEqual(BENCHMARK.MatchResult("CAUSAL", 21, 21), match)
        self.assertFalse(BENCHMARK.hit_at_cutoff(match, 20))

    def test_rankvar_converts_variant_ranks_to_unique_gene_ranks(self):
        path = self.write(
            self.root / "sample.rank_var.tsv",
            "Gene.refGene\tpathogenecity_score\trank\n"
            "GENEB\t0.8\t2.0\n"
            "GENEA\t0.9\t1.0\n"
            "genea\t0.9\t1.0\n"
            "GENEC;GENED\t0.7\t3.0\n",
        )
        rankings = BENCHMARK.read_rankvar_ranking(path)
        self.assertEqual(
            [
                ("GENEA", 1, 1),
                ("GENEB", 2, 2),
                ("GENEC", 3, 3),
                ("GENED", 4, 3),
            ],
            [(r.gene, r.rank, r.raw_rank) for r in rankings],
        )

    def test_rankvar_rejects_non_integral_dense_rank(self):
        path = self.write(
            self.root / "bad.rank_var.tsv",
            "Gene.refGene\trank\nGENEA\t1.5\n",
        )
        with self.assertRaisesRegex(BENCHMARK.RankingError, "non-integral rank"):
            BENCHMARK.read_rankvar_ranking(path)

    def test_cli_uses_norm_alt_pipevar_filename_by_default(self):
        args = BENCHMARK.parse_args(
            [
                "--truth-tsv",
                "truth.tsv",
                "--pipevar-dir",
                "pipevar",
                "--exomiser-dir",
                "exomiser",
                "--output-prefix",
                "benchmark",
            ]
        )
        self.assertEqual("{sample}.norm.alt.prio_gene.vcf", args.pipevar_template)
        self.assertIsNone(args.rankvar_dir)
        self.assertEqual("{sample}.norm.alt.rank_var.tsv", args.rankvar_template)

    def test_best_match_uses_best_of_multiple_truth_genes(self):
        rankings = (
            BENCHMARK.RankedGene("GENE1", 10, 12),
            BENCHMARK.RankedGene("GENE2", 5, 6),
        )
        match = BENCHMARK.find_best_match(("GENE1", "GENE2"), rankings)
        self.assertEqual(BENCHMARK.MatchResult("GENE2", 5, 6), match)
        self.assertFalse(BENCHMARK.hit_at_cutoff(match, 1))
        self.assertTrue(BENCHMARK.hit_at_cutoff(match, 5))
        self.assertTrue(BENCHMARK.hit_at_cutoff(match, 10))

    def test_top20_cutoff_includes_rank_20_but_not_rank_21(self):
        rank_20 = BENCHMARK.MatchResult("GENE20", 20, 20)
        rank_21 = BENCHMARK.MatchResult("GENE21", 21, 21)
        self.assertTrue(BENCHMARK.hit_at_cutoff(rank_20, 20))
        self.assertFalse(BENCHMARK.hit_at_cutoff(rank_21, 20))

    def test_end_to_end_reports_use_paired_complete_denominator(self):
        truth = self.write(
            self.root / "truth.tsv",
            "sample_id\tgenes\nS1\tGENE1\nS2\tGENE2\nS3\tGENE3\n",
        )
        self.write(
            self.pipevar_dir / "S1.norm.alt.prio_gene.vcf",
            "#RANK\tGENE\tVCF_LINE\n1\tGENE1\tchr1\t1\n",
        )
        self.write(
            self.exomiser_dir / "S1.exomiser.input.vcf.gz-results.genes.tsv",
            "#RANK\tGENE_SYMBOL\n1\tOTHER\n2\tGENE1\n",
        )
        self.write(
            self.pipevar_dir / "S1.norm.alt.rank_var.tsv",
            "Gene.refGene\tpathogenecity_score\trank\n"
            "OTHER\t0.9\t1.0\nGENE1\t0.8\t2.0\n",
        )
        self.write(
            self.pipevar_dir / "S2.norm.alt.prio_gene.vcf",
            "#RANK\tGENE\tVCF_LINE\n",
        )
        self.write(
            self.exomiser_dir / "S2.exomiser.input.vcf.gz-results.genes.tsv",
            "#RANK\tGENE_SYMBOL\n1\tGENE2\n",
        )
        self.write(
            self.pipevar_dir / "S2.norm.alt.rank_var.tsv",
            "Gene.refGene\tpathogenecity_score\trank\n"
            "OTHER\t0.9\t1.0\nGENE2\t0.8\t2.0\n",
        )
        self.write(
            self.pipevar_dir / "S3.norm.alt.prio_gene.vcf",
            "#RANK\tWRONG\n1\tGENE3\n",
        )
        output_prefix = self.root / "benchmark"
        args = argparse.Namespace(
            truth_tsv=truth,
            pipevar_dir=self.pipevar_dir,
            rankvar_dir=None,
            exomiser_dir=self.exomiser_dir,
            output_prefix=output_prefix,
            pipevar_template="{sample}.norm.alt.prio_gene.vcf",
            rankvar_template="{sample}.norm.alt.rank_var.tsv",
            exomiser_template="{sample}.exomiser.input.vcf.gz-results.genes.tsv",
        )
        self.assertEqual(0, BENCHMARK.run(args))

        with Path(f"{output_prefix}.summary.tsv").open(newline="") as handle:
            summary = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(["1", "5", "10", "20"], [row["cutoff"] for row in summary])
        self.assertEqual("2", summary[0]["paired_samples"])
        self.assertEqual("1", summary[0]["pipevar_hits"])
        self.assertEqual("0", summary[0]["rankvar_hits"])
        self.assertEqual("1", summary[0]["exomiser_hits"])
        self.assertEqual("1", summary[0]["hit_difference_pipevar_minus_rankvar"])
        self.assertEqual("1", summary[0]["pipevar_only_vs_rankvar"])
        self.assertEqual("1", summary[0]["pipevar_only"])
        self.assertEqual("1", summary[0]["exomiser_only"])

        with Path(f"{output_prefix}.details.tsv").open(newline="") as handle:
            details = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual("empty", details[1]["pipevar_status"])
        self.assertEqual("malformed", details[2]["pipevar_status"])
        self.assertEqual("missing", details[2]["rankvar_status"])
        self.assertEqual("missing", details[2]["exomiser_status"])
        self.assertEqual("0", details[2]["included_in_paired_summary"])

        with Path(f"{output_prefix}.not_found.tsv").open(newline="") as handle:
            not_found = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(["S2", "S3"], [row["sample_id"] for row in not_found])
        self.assertEqual("pipevar", not_found[0]["not_in_top20"])
        self.assertEqual("0", not_found[0]["pipevar_truth_gene_top20"])
        self.assertEqual("1", not_found[0]["exomiser_truth_gene_top20"])
        self.assertEqual(
            "0", not_found[0]["pipevar_truth_gene_found_anywhere"]
        )
        self.assertEqual(
            "1", not_found[0]["exomiser_truth_gene_found_anywhere"]
        )
        self.assertEqual("GENE2", not_found[0]["exomiser_matched_gene"])
        self.assertEqual("1", not_found[0]["exomiser_best_rank"])
        self.assertEqual("1:GENE2", not_found[0]["exomiser_top20_genes"])
        self.assertEqual("both", not_found[1]["not_in_top20"])
        self.assertEqual("malformed", not_found[1]["pipevar_status"])
        self.assertEqual("missing", not_found[1]["exomiser_status"])

    def test_not_found_report_excludes_samples_found_by_both_tools(self):
        truth = BENCHMARK.TruthSample("S1", ("GENE1",))
        ranking = (BENCHMARK.RankedGene("GENE1", 1, 1),)
        valid = BENCHMARK.RankingResult("ok", ranking)
        result = BENCHMARK.SampleResult(
            truth=truth,
            pipevar_path=Path("pipevar"),
            rankvar_path=Path("rankvar"),
            exomiser_path=Path("exomiser"),
            pipevar_result=valid,
            rankvar_result=valid,
            exomiser_result=valid,
            pipevar_match=BENCHMARK.MatchResult("GENE1", 1, 1),
            rankvar_match=BENCHMARK.MatchResult("GENE1", 1, 1),
            exomiser_match=BENCHMARK.MatchResult("GENE1", 1, 1),
        )
        report = self.root / "not_found.tsv"
        self.assertEqual(0, BENCHMARK.write_not_found_report(report, [result]))
        with report.open(newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual([], rows)

    def test_not_found_category_distinguishes_each_tool_and_both(self):
        truth = BENCHMARK.TruthSample("S1", ("GENE1",))
        valid = BENCHMARK.RankingResult("ok", ())

        def make(pipevar_rank, exomiser_rank):
            return BENCHMARK.SampleResult(
                truth, Path(), Path(), Path(), valid, valid, valid,
                BENCHMARK.MatchResult(rank=pipevar_rank),
                BENCHMARK.MatchResult(),
                BENCHMARK.MatchResult(rank=exomiser_rank),
            )

        self.assertEqual("pipevar", BENCHMARK.not_found_category(make(None, 2)))
        self.assertEqual("exomiser", BENCHMARK.not_found_category(make(2, None)))
        self.assertEqual("both", BENCHMARK.not_found_category(make(None, None)))
        self.assertEqual("", BENCHMARK.not_found_category(make(2, 2)))

    def test_not_found_report_treats_rank_21_as_a_top20_miss(self):
        truth = BENCHMARK.TruthSample("S1", ("SETBP1",))
        pipevar = BENCHMARK.RankingResult("ok", ())
        exomiser_ranking = (BENCHMARK.RankedGene("SETBP1", 21, 21),)
        exomiser = BENCHMARK.RankingResult("ok", exomiser_ranking)
        result = BENCHMARK.SampleResult(
            truth=truth,
            pipevar_path=Path("pipevar"),
            rankvar_path=Path("rankvar"),
            exomiser_path=Path("exomiser"),
            pipevar_result=pipevar,
            rankvar_result=pipevar,
            exomiser_result=exomiser,
            pipevar_match=BENCHMARK.MatchResult(),
            rankvar_match=BENCHMARK.MatchResult(),
            exomiser_match=BENCHMARK.MatchResult("SETBP1", 21, 21),
        )

        self.assertEqual("both", BENCHMARK.not_found_category(result))
        report = self.root / "rank21.not_found.tsv"
        self.assertEqual(1, BENCHMARK.write_not_found_report(report, [result]))
        with report.open(newline="") as handle:
            row = next(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual("both", row["not_in_top20"])
        self.assertEqual("0", row["exomiser_truth_gene_top20"])
        self.assertEqual("1", row["exomiser_truth_gene_found_anywhere"])
        self.assertEqual("SETBP1", row["exomiser_matched_gene"])
        self.assertEqual("21", row["exomiser_best_rank"])
        self.assertNotIn("SETBP1", row["exomiser_top20_genes"])

    def test_overridden_templates_resolve_nested_paths(self):
        sample = BENCHMARK.TruthSample("S1", ("GENE1",))
        (self.pipevar_dir / "S1").mkdir()
        (self.exomiser_dir / "S1").mkdir()
        self.write(
            self.pipevar_dir / "S1" / "rank.tsv", "#RANK\tGENE\n1\tGENE1\n"
        )
        self.write(
            self.exomiser_dir / "S1" / "genes.tsv",
            "#RANK\tGENE_SYMBOL\n1\tGENE1\n",
        )
        self.write(
            self.pipevar_dir / "S1" / "rankvar.tsv",
            "Gene.refGene\trank\nGENE1\t1.0\n",
        )
        results = BENCHMARK.evaluate_samples(
            [sample],
            self.pipevar_dir,
            self.exomiser_dir,
            "{sample}/rank.tsv",
            "{sample}/genes.tsv",
            self.pipevar_dir,
            "{sample}/rankvar.tsv",
        )
        self.assertTrue(results[0].paired)

    def test_zero_paired_samples_writes_reports_and_returns_one(self):
        truth = self.write(
            self.root / "truth.tsv", "sample_id\tgenes\nMISSING\tGENE1\n"
        )
        args = argparse.Namespace(
            truth_tsv=truth,
            pipevar_dir=self.pipevar_dir,
            rankvar_dir=None,
            exomiser_dir=self.exomiser_dir,
            output_prefix=self.root / "none",
            pipevar_template="{sample}.norm.alt.prio_gene.vcf",
            rankvar_template="{sample}.norm.alt.rank_var.tsv",
            exomiser_template="{sample}.exomiser.input.vcf.gz-results.genes.tsv",
        )
        self.assertEqual(1, BENCHMARK.run(args))
        with (self.root / "none.summary.tsv").open(newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertTrue(all(row["pipevar_rate"] == "NA" for row in rows))
        self.assertTrue((self.root / "none.not_found.tsv").is_file())


if __name__ == "__main__":
    unittest.main()
