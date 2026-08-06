import re
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
DOCKER_WORK = REPO.parent / "docker_work"


class RepairContractTests(unittest.TestCase):
    def test_phenosv_modules_partition_merge_and_warn_in_light_mode(self):
        for relative in ("modules/phenosv/main.nf", "modules/multi_phenosv/main.nf"):
            source = (REPO / relative).read_text(encoding="utf-8")
            self.assertIn("canonical_bed", source)
            self.assertIn("phenosv_bedpe", source)
            self.assertIn("members_tsv", source)
            self.assertIn("--members $members_tsv", source)
            self.assertIn("PhenoSV-light has reduced accuracy for translocations", source)
            self.assertEqual(1, source.count("head -n 1 ${out_prefix}.phenosv.simple.tsv"))
            self.assertIn("tail -n +2 ${out_prefix}.phenosv.bnd.tsv", source)

    def test_phenosv_image_does_not_patch_upstream_checkout(self):
        dockerfile = (DOCKER_WORK / "phenosv/Dockerfile").read_text(encoding="utf-8")
        self.assertIn("phenosv_0.2", dockerfile)
        self.assertNotRegex(dockerfile, r"sed .*PhenoSV|operation_function|git apply|patch ")

    def test_rankscore_ad_is_wired_to_all_callers(self):
        calls = []
        for path in (REPO / "subworkflows").glob("*/main.nf"):
            source = path.read_text(encoding="utf-8")
            calls.extend(re.findall(r"(?:Rankscore_analysis|multi_rankscore)\([^\n]+\)", source))
        self.assertEqual(25, len(calls))
        self.assertTrue(all(re.search(r",\s*ad\s*,\s*phen2gene_top_n\)$", call) for call in calls))

        for relative in ("modules/rankscore_analysis/main.nf", "modules/multi_rankscore/main.nf"):
            source = (REPO / relative).read_text(encoding="utf-8")
            self.assertIn("rankscore_0.2.22", source)
            self.assertRegex(source, r"val ad[\s\S]+?\$gq \$ad \$rankscore_softwares")

    def test_reserved_image_pins_are_scoped_to_consumers(self):
        prioritization = {
            "longphase", "multi_longphase", "ngs_prio", "multi_ngs_prio",
            "snp_prio", "multi_snp_prio", "sv_prio", "multi_sv_prio",
        }
        pinned = {
            path.parent.name
            for path in (REPO / "modules").glob("*/main.nf")
            if "longphase_0.2.33" in path.read_text(encoding="utf-8")
        }
        self.assertEqual(prioritization, pinned)
        for merge_module in (
            "merge_longread_sv_callers", "multi_merge_longread_sv_callers",
            "merge_shortread_sv_callers", "multi_merge_shortread_sv_callers",
        ):
            self.assertNotIn(
                "longphase_0.2.33",
                (REPO / "modules" / merge_module / "main.nf").read_text(encoding="utf-8"),
            )

    def test_rankscore_clinvar_filter_is_header_indexed(self):
        source = (DOCKER_WORK / "rankscore/clinvar.sh").read_text(encoding="utf-8")
        self.assertIn('h["CLNSIG"]', source)
        self.assertIn("accepted_clnsig", source)
        self.assertNotIn("&& /Pathogenic/", source)


if __name__ == "__main__":
    unittest.main()
