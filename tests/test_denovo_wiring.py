import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


class DenovoWiringTests(unittest.TestCase):
    called_snv_workflows = [
        "subworkflows/input_csv_alignment_vcf_snp/main.nf",
        "subworkflows/input_csv_alignment_ngs_snp/main.nf",
        "subworkflows/input_csv_alignment_long_snp/main.nf",
        "subworkflows/input_csv_alignment_all_ngs/main.nf",
        "subworkflows/input_csv_alignment_all_longphase/main.nf",
    ]
    annotated_snv_workflows = [
        "subworkflows/annotated_snv_prio_core/main.nf",
        "subworkflows/input_csv_annotated_all_ngs/main.nf",
        "subworkflows/input_csv_annotated_snv_called_sv_ngs/main.nf",
    ]
    sv_workflows = [
        "subworkflows/input_csv_alignment_vcf_sv/main.nf",
        "subworkflows/input_csv_alignment_ngs_sv/main.nf",
        "subworkflows/input_csv_alignment_long_sv/main.nf",
        "subworkflows/input_csv_alignment_all_ngs/main.nf",
        "subworkflows/input_csv_alignment_all_longphase/main.nf",
        "subworkflows/input_csv_annotated_all_ngs/main.nf",
        "subworkflows/input_csv_annotated_snv_called_sv_ngs/main.nf",
    ]

    def read(self, relative):
        return (ROOT / relative).read_text()

    def test_called_snv_workflows_filter_raw_vcf_before_annovar(self):
        for relative in self.called_snv_workflows:
            with self.subTest(relative=relative):
                text = self.read(relative)
                self.assertIn("include { DENOVO_SNV_VCF_FILTER_CORE }", text)
                self.assertLess(text.index("DENOVO_SNV_VCF_FILTER_CORE("), text.index("multi_annovar(annovar_input)"))
                self.assertTrue("snp_for_annotation" in text or "snv_for_annotation" in text)

    def test_imported_annovar_workflows_keep_paired_txt_vcf_filtering(self):
        for relative in self.annotated_snv_workflows:
            with self.subTest(relative=relative):
                text = self.read(relative)
                self.assertIn("include { DENOVO_SNV_FILTER_CORE }", text)
                self.assertIn("validated_for_downstream", text)

    def test_all_active_sv_workflows_filter_before_common_sv_filter(self):
        for relative in self.sv_workflows:
            with self.subTest(relative=relative):
                text = self.read(relative)
                denovo_at = text.index("DENOVO_SV_FILTER_CORE(")
                annovar_at = min(
                    pos for pos in (
                        text.find("multi_annovar_sv(sv_result_annovar)"),
                        text.find("multi_annovar_sv(sniffles_result_annovar)"),
                    ) if pos >= 0
                )
                common_at = text.index("multi_common_sv_filter(annovar_sv_for_downstream)")
                self.assertLess(denovo_at, annovar_at)
                self.assertLess(denovo_at, common_at)
                self.assertIn("sv_for_annotation", text)

    def test_top_level_validates_and_passes_pedigree_settings(self):
        text = self.read("main.nf")
        self.assertIn("def clean_denovo_filter", text)
        self.assertIn("family '${familyId}' is missing a proband", text)
        self.assertIn("denovo_pedigree = Channel.fromPath(params.input_csv).first()", text)
        self.assertGreaterEqual(text.count("clean_denovo_filter, denovo_pedigree"), 11)
        self.assertIn("def parseCsvRecord", text)
        self.assertIn("duplicate sample '${sample}' in de novo CSV", text)
        self.assertIn("role != 'proband'", text)
        self.assertIn("requested --sex_column", text)

    def test_adapters_collect_complete_keyed_records(self):
        snv = self.read("subworkflows/denovo_snv_filter_core/main.nf")
        sv = self.read("subworkflows/denovo_sv_filter_core/main.nf")
        self.assertIn("annovar_records.collect()", snv)
        self.assertIn("ordered.collect { it[0] }", snv)
        self.assertIn("ordered.collect { it[1] }", snv)
        self.assertIn("ordered.collect { it[2] }", snv)
        self.assertIn("sv_records.collect()", sv)
        self.assertIn("ordered.collect { it[0] }", sv)
        self.assertIn("ordered.collect { it[1] }", sv)
        self.assertIn("filtered.bindings", snv)
        self.assertIn("filtered.bindings", sv)
        self.assertNotIn("replaceFirst", snv)
        self.assertNotIn("replaceFirst", sv)

    def test_both_modules_use_the_same_image_contract(self):
        snv = self.read("modules/multi_denovo_snv_filter/main.nf")
        raw_snv = self.read("modules/multi_denovo_snv_vcf_filter/main.nf")
        sv = self.read("modules/multi_denovo_sv_filter/main.nf")
        tag = "beoungl/docker_test:denovo_snv_sv_filter_0.3"
        self.assertIn(tag, snv)
        self.assertIn(tag, raw_snv)
        self.assertIn(tag, sv)
        self.assertIn("--snv-input", snv)
        self.assertIn("--snv-vcf-input", raw_snv)
        self.assertIn("--sv-input", sv)

    def test_proband_sex_is_rejoined_after_family_filtering(self):
        workflows = self.called_snv_workflows + self.sv_workflows + self.annotated_snv_workflows
        for relative in workflows:
            with self.subTest(relative=relative):
                text = self.read(relative)
                self.assertIn("proband_keys", text)
                self.assertIn("input_meta", text)
                self.assertIn("sex", text)
                self.assertIn("failOnMismatch: true", text)

    def test_publish_pattern_covers_outputs_and_summaries(self):
        config = self.read("nextflow.config")
        self.assertIn("multi_denovo_snv_vcf_filter", config)
        self.assertIn("denovo.*.bindings.tsv", config)


if __name__ == "__main__":
    unittest.main()
