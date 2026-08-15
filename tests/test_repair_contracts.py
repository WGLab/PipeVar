import re
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
DOCKER_WORK = REPO.parent / "docker_work"


class RepairContractTests(unittest.TestCase):
    def test_survivor_single_and_batch_use_embedded_phenosv_converter(self):
        dockerfile = (DOCKER_WORK / "survivor/Dockerfile").read_text(encoding="utf-8")
        self.assertIn("COPY prepare_phenosv_inputs.py /prepare_phenosv_inputs.py", dockerfile)
        for relative in ("modules/survivor/main.nf", "modules/multi_survivor/main.nf"):
            source = (REPO / relative).read_text(encoding="utf-8")
            self.assertIn("beoungl/docker_test:survivor_0.2", source)
            self.assertIn("python3 /prepare_phenosv_inputs.py", source)

    def test_phenosv_modules_partition_merge_and_warn_in_light_mode(self):
        for relative in ("modules/phenosv/main.nf", "modules/multi_phenosv/main.nf"):
            source = (REPO / relative).read_text(encoding="utf-8")
            self.assertIn("canonical_bed", source)
            self.assertIn("phenosv_bedpe", source)
            self.assertIn("members_tsv", source)
            self.assertIn("--members $members_tsv", source)
            self.assertIn("PhenoSV-light has reduced accuracy for translocations", source)
            self.assertIn("PHENOSV_EVENT_ID", source)
            self.assertIn("PHENOSV_GENE_SCORE", source)
            self.assertEqual(1, source.count("head -n 1 ${out_prefix}.phenosv.simple.tsv"))
            self.assertIn("tail -n +2 ${out_prefix}.phenosv.bnd.tsv", source)

    def test_phenosv_image_does_not_patch_upstream_checkout(self):
        dockerfile = (DOCKER_WORK / "phenosv/Dockerfile").read_text(encoding="utf-8")
        self.assertIn("phenosv_0.3", dockerfile)
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
            self.assertIn("rankscore_0.3.0", source)
            self.assertRegex(source, r"val ad[\s\S]+?\$gq \$ad \$rankscore_softwares")

    def test_reserved_image_pins_are_scoped_to_consumers(self):
        prioritization = {
            "longphase", "multi_longphase", "ngs_prio", "multi_ngs_prio",
            "snp_prio", "multi_snp_prio", "sv_prio", "multi_sv_prio",
        }
        pinned = {
            path.parent.name
            for path in (REPO / "modules").glob("*/main.nf")
            if "longphase_0.4.0" in path.read_text(encoding="utf-8")
        }
        self.assertEqual(prioritization, pinned)
        for merge_module in (
            "merge_longread_sv_callers", "multi_merge_longread_sv_callers",
            "merge_shortread_sv_callers", "multi_merge_shortread_sv_callers",
        ):
            self.assertNotIn(
                "longphase_0.4.0",
                (REPO / "modules" / merge_module / "main.nf").read_text(encoding="utf-8"),
            )

    def test_common_sv_filter_uses_canonical_duplication_image(self):
        for module in ("common_sv_filter", "multi_common_sv_filter"):
            source = (REPO / "modules" / module / "main.nf").read_text(encoding="utf-8")
            self.assertIn("common_sv_filter_0.5", source)
        dockerfile = (DOCKER_WORK / "common_sv_filter/Dockerfile").read_text(encoding="utf-8")
        self.assertIn("common_sv_filter_0.5", dockerfile)

    def test_rankscore_clinvar_filter_is_header_indexed(self):
        source = (DOCKER_WORK / "rankscore/clinvar.sh").read_text(encoding="utf-8")
        self.assertIn('h["CLNSIG"]', source)
        self.assertIn("accepted_clnsig", source)
        self.assertNotIn("&& /Pathogenic/", source)

    def test_preannotated_validator_is_copy_free_and_requires_gt_but_not_ps(self):
        module = (REPO / "modules/validate_preannotated_annovar_pair/main.nf").read_text(encoding="utf-8")
        validator = (
            DOCKER_WORK
            / "validate_preannotated_annovar_pair/validate_preannotated_annovar_pair.py"
        ).read_text(encoding="utf-8")
        dockerfile = (
            DOCKER_WORK / "validate_preannotated_annovar_pair/Dockerfile"
        ).read_text(encoding="utf-8")

        self.assertIn("validate_preannotated_annovar_pair_0.4", module)
        self.assertIn("validate_preannotated_annovar_pair_0.4", dockerfile)
        self.assertIn("procps", dockerfile)
        self.assertIn('if "GT" not in format_fields', validator)
        self.assertIn("PS is optional", validator)
        self.assertNotIn("shutil", validator)
        self.assertNotIn("copy_file", validator)
        self.assertNotIn("--validated-txt", validator)
        self.assertNotIn("--validated-vcf", validator)
        self.assertIn("process validate_preannotated_annovar_pair", module)
        self.assertIn("tuple val(out_prefix), path(annovar_txt), path(annovar_vcf)", module)
        self.assertNotIn("validate_preannotated_annovar_pair_check", module)
        self.assertNotIn("validation_ok", module)
        self.assertNotIn(".join(", module)
        self.assertIn("stub:", module)
        self.assertNotIn(".validated.hg38_multianno", module)
        config = (REPO / "nextflow.config").read_text(encoding="utf-8")
        self.assertIn("withName: 'validate_preannotated_annovar_pair'", config)
        self.assertNotIn("withName: 'validate_preannotated_annovar_pair_check'", config)

if __name__ == "__main__":
    unittest.main()
