import re
import unittest
from pathlib import Path


PIPEVAR = Path(__file__).resolve().parents[1]


def read(relative_path):
    return (PIPEVAR / relative_path).read_text()


class LightCallerWiringTests(unittest.TestCase):
    def test_light_is_normalized_and_validated_strictly(self):
        main = read("main.nf")

        self.assertIn(
            "def clean_light  = params.light == null ? 'no' : "
            "params.light.toString().trim().toLowerCase()",
            main,
        )
        self.assertIn("if (!valid_yes_no.contains(clean_light))", main)
        self.assertIn('You provided: --light "${params.light}"', main)
        self.assertIn("--light yes", main)
        self.assertIn("--light no", main)

    def test_caller_selectors_are_plain_scalars(self):
        main = read("main.nf")

        self.assertIn("def input_bam = null", main)
        self.assertIn(
            "def short_snp_caller = clean_light == 'yes' ? "
            "'haplotypecaller' : 'deepvariant'",
            main,
        )
        self.assertIn(
            "def long_snp_caller = clean_light == 'yes' ? "
            "'nanocaller' : 'clair3'",
            main,
        )
        self.assertNotRegex(main, r"short_snp_caller\s*=\s*Channel\.value")
        self.assertNotRegex(main, r"long_snp_caller\s*=\s*Channel\.value")

    def test_short_read_references_are_separate(self):
        main = read("main.nf")

        canonical_reference = re.compile(
            r"ref_fa\s*=\s*Channel.*?"
            r"return \[ fa_file, fai_file \].*?"
            r"\.first\(\)",
            re.DOTALL,
        )
        gatk_reference = re.compile(
            r"gatk_ref_fa\s*=\s*Channel.*?"
            r"def dict_file = .*?"
            r"return \[ fa_file, fai_file, dict_file \].*?"
            r"\.first\(\)",
            re.DOTALL,
        )
        self.assertRegex(main, canonical_reference)
        self.assertRegex(main, gatk_reference)

        for workflow_call in (
            "INPUT_CSV_NGS_SNP(input_bam, input_meta, ref_fa, gatk_ref_fa,",
            "INPUT_CSV_ALIGNMENT_ALL_NGS(input_bam, input_meta, ref_fa, gatk_ref_fa,",
            "SINGLE_ALIGNMENT_NGS_SNP(bam, out_prefix, ref_fa, gatk_ref_fa,",
            "SINGLE_ALIGNMENT_ALL_NGS(bam, out_prefix, ref_fa, gatk_ref_fa,",
        ):
            self.assertIn(workflow_call, main)

    def test_all_consolidated_workflows_branch_on_resolved_caller_name(self):
        expected = {
            "subworkflows/single_alignment_ngs_snp/main.nf": (
                'caller_mode == "haplotypecaller"',
                "haplotypecaller.out : deepvariant.out",
            ),
            "subworkflows/single_alignment_all_ngs/main.nf": (
                'caller_mode == "haplotypecaller"',
                "haplotypecaller.out : deepvariant.out",
            ),
            "subworkflows/input_csv_alignment_ngs_snp/main.nf": (
                'caller_mode == "haplotypecaller"',
                "multi_haplotypecaller(",
            ),
            "subworkflows/input_csv_alignment_all_ngs/main.nf": (
                'caller_mode == "haplotypecaller"',
                "multi_haplotypecaller(",
            ),
            "subworkflows/single_alignment_long_snp/main.nf": (
                'caller_mode == "nanocaller"',
                "nanocaller.out : clair3.out",
            ),
            "subworkflows/single_alignment_all_longphase/main.nf": (
                'caller_mode == "nanocaller"',
                "nanocaller.out : clair3.out",
            ),
            "subworkflows/input_csv_alignment_long_snp/main.nf": (
                'caller_mode == "nanocaller"',
                "multi_nanocaller(",
            ),
            "subworkflows/input_csv_alignment_all_longphase/main.nf": (
                'caller_mode == "nanocaller"',
                "multi_nanocaller(",
            ),
        }

        for relative_path, required_fragments in expected.items():
            workflow = read(relative_path)
            for fragment in required_fragments:
                self.assertIn(fragment, workflow, relative_path)

    def test_haplotypecaller_uses_only_gatk_reference(self):
        for relative_path in (
            "subworkflows/single_alignment_ngs_snp/main.nf",
            "subworkflows/single_alignment_all_ngs/main.nf",
            "subworkflows/input_csv_alignment_ngs_snp/main.nf",
            "subworkflows/input_csv_alignment_all_ngs/main.nf",
        ):
            workflow = read(relative_path)
            self.assertIn("gatk_ref_fa", workflow, relative_path)
            self.assertNotRegex(
                workflow,
                r"(?:multi_)?haplotypecaller\([^\n]*,\s*ref_fa(?:,|\))",
                relative_path,
            )

    def test_batch_regions_are_keyed_and_staged_with_caller_input(self):
        workflow_paths = (
            "subworkflows/input_csv_alignment_ngs_snp/main.nf",
            "subworkflows/input_csv_alignment_all_ngs/main.nf",
            "subworkflows/input_csv_alignment_long_snp/main.nf",
            "subworkflows/input_csv_alignment_all_longphase/main.nf",
        )
        for relative_path in workflow_paths:
            workflow = read(relative_path)
            self.assertIn("tuple(out_prefix, [])", workflow, relative_path)
            self.assertIn(
                ".join(caller_regions, failOnMismatch: true, failOnDuplicate: true)",
                workflow,
                relative_path,
            )

        module_paths = (
            "modules/multi_deepvariant/main.nf",
            "modules/multi_haplotypecaller/main.nf",
            "modules/multi_clair3/main.nf",
            "modules/multi_nanocaller/main.nf",
        )
        for relative_path in module_paths:
            module = read(relative_path)
            self.assertRegex(
                module,
                r"tuple val\(out_prefix\),[^\n]*path\(bed_file\)",
                relative_path,
            )
            self.assertNotRegex(module, r"(?m)^\s*val bed_file\s*$", relative_path)

    def test_keyed_join_is_order_independent(self):
        alignments = [
            ("sample_a", "a.bam", "a.bam.bai"),
            ("sample_b", "b.bam", "b.bam.bai"),
        ]
        beds_reversed = [
            ("sample_b", "b.bed"),
            ("sample_a", "a.bed"),
        ]
        beds_by_sample = dict(beds_reversed)

        joined = [
            sample_tuple + (beds_by_sample[sample_tuple[0]],)
            for sample_tuple in alignments
        ]

        self.assertEqual(
            joined,
            [
                ("sample_a", "a.bam", "a.bam.bai", "a.bed"),
                ("sample_b", "b.bam", "b.bam.bai", "b.bed"),
            ],
        )

    def test_definite_caller_module_defects_are_fixed(self):
        haplotypecaller = read("modules/multi_haplotypecaller/main.nf")
        deepvariant = read("modules/deepvariant/main.nf")

        self.assertIn(
            "--resource:dbsnp,known=true,training=false,truth=false,prior=2.0 "
            "$dbsnp_vcf",
            haplotypecaller,
        )
        self.assertNotIn(
            "--resource:dbsnp,known=true,training=false,truth=false,prior=2.0 "
            "$hapmap_vcf",
            haplotypecaller,
        )
        self.assertIn("tuple path(bam), path(bam_index)", deepvariant)
        self.assertIn("tuple path(ref_fa), path(fa_index)", deepvariant)
        self.assertNotRegex(deepvariant, r"tuple path\([^)]*\), path\(index\)")


if __name__ == "__main__":
    unittest.main()
