import csv
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
GENERATOR = ROOT / "scripts" / "generate_input_csv.sh"


class AnnotatedSnvSvNoAlignmentTests(unittest.TestCase):
    def test_generator_builds_snv_sv_manifest_without_alignment(self):
        with tempfile.TemporaryDirectory() as tmp_name:
            tmp = Path(tmp_name)
            phenotype_dir = tmp / "phenotype"
            snv_txt_dir = tmp / "snv_txt"
            snv_vcf_dir = tmp / "snv_vcf"
            sv_vcf_dir = tmp / "sv_vcf"
            for directory in (phenotype_dir, snv_txt_dir, snv_vcf_dir, sv_vcf_dir):
                directory.mkdir()
            (phenotype_dir / "case1.hpo.txt").write_text("HP:0001250\n", encoding="utf-8")
            (snv_txt_dir / "case1.hg38_multianno.txt").write_text("Chr\tStart\n", encoding="utf-8")
            (snv_vcf_dir / "case1.hg38_multianno.vcf").write_text("##fileformat=VCFv4.2\n", encoding="utf-8")
            expected_sv = sv_vcf_dir / "case1.hg38_multianno.vcf"
            expected_sv.write_text("##fileformat=VCFv4.2\n", encoding="utf-8")
            output_csv = tmp / "generated.csv"

            answers = "\n".join([
                "2", "1", "2", str(phenotype_dir), ".hpo.txt", str(output_csv),
                "n", "n", str(snv_txt_dir), "", str(snv_vcf_dir), "", "y", "n",
                str(sv_vcf_dir), "", "",
            ])
            result = subprocess.run(
                ["bash", str(GENERATOR)],
                input=answers,
                text=True,
                capture_output=True,
                check=False,
            )
            self.assertEqual(0, result.returncode, result.stderr)
            self.assertIn("Matched rows        : 1", result.stdout)
            with output_csv.open(newline="", encoding="utf-8") as handle:
                row = next(csv.DictReader(handle))
            self.assertEqual(str(expected_sv), row["sv_vcf_path"])
            self.assertEqual("", row["alignment_path"])
            self.assertEqual("", row["alignment_index_path"])

    def test_updater_accepts_blank_alignment_and_reports_skip_reasons(self):
        with tempfile.TemporaryDirectory() as tmp_name:
            tmp = Path(tmp_name)
            sv_dir = tmp / "sv"
            sv_dir.mkdir()
            matched_sv = sv_dir / "case1.hg38_multianno.vcf"
            matched_sv.write_text("##fileformat=VCFv4.2\n", encoding="utf-8")

            input_csv = tmp / "input.csv"
            output_csv = tmp / "updated.csv"
            header = [
                "sample", "input_kind", "phenotype_path", "phenotype_format",
                "age_of_onset", "snv_txt_path", "snv_vcf_path", "sv_vcf_path",
                "vcf_path", "alignment_path", "alignment_index_path",
            ]
            rows = [
                ["case1", "ANNOTATED_SNV", "/note/1", "hpo", "", "/snv/1.txt", "/snv/1.vcf", "", "", "", ""],
                ["case2", "annotated_snv", "/note/2", "hpo", "", "/snv/2.txt", "/snv/2.vcf", "/existing/2.vcf", "", "", ""],
                ["case3", "vcf_snv", "/note/3", "hpo", "", "", "", "", "/snv/3.vcf", "", ""],
                ["case4", "annotated_snv", "/note/4", "hpo", "", "/snv/4.txt", "/snv/4.vcf", "", "", "", ""],
            ]
            with input_csv.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.writer(handle, lineterminator="\n")
                writer.writerow(header)
                writer.writerows(rows)

            answers = "\n".join([
                "3", str(input_csv), str(output_csv), str(sv_dir), "", "n", "",
            ])
            result = subprocess.run(
                ["bash", str(GENERATOR)],
                input=answers,
                text=True,
                capture_output=True,
                check=False,
            )
            self.assertEqual(0, result.returncode, result.stderr)
            self.assertIn("SV paths updated : 1", result.stdout)
            self.assertIn("Non-annotated    : 1", result.stdout)
            self.assertIn("Existing kept    : 1", result.stdout)
            self.assertIn("Missing matches  : 1", result.stdout)

            with output_csv.open(newline="", encoding="utf-8") as handle:
                updated = list(csv.DictReader(handle))
            self.assertEqual(str(matched_sv), updated[0]["sv_vcf_path"])
            self.assertEqual("", updated[0]["alignment_path"])
            self.assertEqual("/existing/2.vcf", updated[1]["sv_vcf_path"])
            self.assertEqual("", updated[3]["sv_vcf_path"])

    def test_updater_rejects_in_place_output(self):
        with tempfile.TemporaryDirectory() as tmp_name:
            input_csv = Path(tmp_name) / "input.csv"
            input_csv.write_text("sample,input_kind\ncase1,annotated_snv\n", encoding="utf-8")
            answers = "\n".join(["3", str(input_csv), str(input_csv), ""])
            result = subprocess.run(
                ["bash", str(GENERATOR)],
                input=answers,
                text=True,
                capture_output=True,
                check=False,
            )
            self.assertNotEqual(0, result.returncode)
            self.assertIn("input and output CSV paths must be different", result.stderr)
            self.assertIn("case1,annotated_snv", input_csv.read_text(encoding="utf-8"))

    def test_pipeline_has_four_homogeneous_annotated_routes(self):
        main = (ROOT / "main.nf").read_text(encoding="utf-8")
        self.assertIn("manifestAnnotatedRowsWithPreannotatedSv > 0 && manifestAnnotatedRowsWithAlignment > 0", main)
        self.assertIn("else if (manifestAnnotatedRowsWithPreannotatedSv > 0)", main)
        self.assertIn("else if (manifestAnnotatedRowsWithAlignment > 0)", main)
        self.assertIn("csv_manifest_mode = 'annotated_snv'", main)
        self.assertIn("input_annotated_snv_sv = Channel", main)
        self.assertIn("INPUT_CSV_ANNOTATED_SNV_SV(", main)
        self.assertIn("SINGLE_ANNOTATED_SNV_SV(", main)
        self.assertIn("must either provide sv_vcf_path for every row or leave it blank", main)
        self.assertIn("--mito yes requires alignment_path", main)
        self.assertNotIn("provides sv_vcf_path without alignment_path", main)

    def test_no_alignment_workflows_exclude_alignment_only_analysis(self):
        batch = (ROOT / "subworkflows/input_csv_annotated_snv_sv/main.nf").read_text(encoding="utf-8")
        single = (ROOT / "subworkflows/single_annotated_snv_sv/main.nf").read_text(encoding="utf-8")
        for source in (batch, single):
            self.assertNotIn("bam_file", source)
            self.assertNotIn("bai_file", source)
            self.assertNotIn("ExpansionHunter", source)
            self.assertNotIn("mito_", source)
            self.assertNotIn("ref_fa", source)
            self.assertIn("ngs_prio", source)
            self.assertIn("variant_html_report_no_repeat", source)
        self.assertIn("DENOVO_SNV_FILTER_CORE", batch)
        self.assertIn("DENOVO_SV_FILTER_CORE", batch)

    def test_generator_prompts_for_sv_and_alignment_independently(self):
        source = GENERATOR.read_text(encoding="utf-8")
        self.assertIn("Include annotated SV VCF paths?", source)
        self.assertIn("Also include BAM/CRAM paths for STR/mito or SV calling?", source)
        self.assertIn('if [[ "$include_sv" == "yes" ]]', source)
        self.assertIn('if [[ "$include_alignment" == "yes" ]]', source)


if __name__ == "__main__":
    unittest.main()
