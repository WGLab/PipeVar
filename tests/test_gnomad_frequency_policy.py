import csv
import importlib.util
import tempfile
import unittest
from pathlib import Path


RANKSCORE_SCRIPT = Path("/home/beoungle/docker_work/rankscore/normalize_gnomad_frequency.py")
RANKVAR_SCRIPT = Path("/home/beoungle/docker_work/rankvar/normalize_gnomad_frequency.py")
PIPEVAR_ROOT = Path(__file__).resolve().parents[1]


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class GnomadFrequencyPolicyTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.rankscore = load_module("rankscore_gnomad_frequency", RANKSCORE_SCRIPT)
        cls.rankvar = load_module("rankvar_gnomad_frequency", RANKVAR_SCRIPT)

    def test_rankscore_and_rankvar_select_identically(self):
        rows = [
            {
                "gnomad41_exome_fafmax_faf95_max": "0.004",
                "gnomad41_genome_fafmax_faf95_max": "0.006",
                "gnomad41_exome_AF_grpmax": "0.02",
                "gnomad41_genome_AF_grpmax": "0.03",
            },
            {
                "gnomad41_exome_fafmax_faf95_max": ".",
                "gnomad41_genome_fafmax_faf95_max": "",
                "gnomad41_exome_AF_grpmax": "9e-3",
                "gnomad41_genome_AF_grpmax": "0.01",
            },
            {},
        ]
        for line_number, row in enumerate(rows, start=2):
            self.assertEqual(
                self.rankscore.select_frequency(row, line_number),
                self.rankvar.select_frequency(row, line_number),
            )
        self.assertEqual((0.006, "GNOMAD41_FAF95_EXOME_GENOME"), self.rankscore.select_frequency(rows[0], 2))
        self.assertEqual((0.01, "GNOMAD41_AF_GRPMAX_EXOME_GENOME"), self.rankscore.select_frequency(rows[1], 3))
        self.assertEqual((0.0, "GNOMAD41_MISSING_AS_ZERO"), self.rankscore.select_frequency(rows[2], 4))

    def test_normalizer_keeps_exact_boundary_and_filters_above(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            source = root / "input.tsv"
            output = root / "output.tsv"
            with source.open("w", newline="") as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=["Chr", "Start", "Ref", "Alt", "gnomad41_exome_AF_grpmax"],
                    delimiter="\t",
                    lineterminator="\n",
                )
                writer.writeheader()
                writer.writerow({"Chr": "1", "Start": "1", "Ref": "A", "Alt": "G", "gnomad41_exome_AF_grpmax": "0.01"})
                writer.writerow({"Chr": "1", "Start": "2", "Ref": "A", "Alt": "T", "gnomad41_exome_AF_grpmax": "0.0101"})
            self.rankscore.normalize_table(source, output, max_af=0.01)
            with output.open(newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(["1"], [row["Start"] for row in rows])

    def test_malformed_nonmissing_frequency_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "malformed"):
            self.rankscore.select_frequency({"gnomad41_exome_AF_grpmax": "not-a-number"}, 9)

    def test_recognized_missing_values_and_absent_columns_normalize_to_zero(self):
        for value in (".", "", "NA", "N/A", "nan", "None"):
            row = {"gnomad41_exome_fafmax_faf95_max": value}
            self.assertEqual(
                (0.0, "GNOMAD41_MISSING_AS_ZERO"),
                self.rankscore.select_frequency(row, 2),
            )
        self.assertEqual(
            (0.0, "GNOMAD41_MISSING_AS_ZERO"),
            self.rankscore.select_frequency({}, 2),
        )

    def test_multiallelic_rows_remain_allele_specific(self):
        rows = [
            {"Alt": "C", "gnomad41_exome_AF_grpmax": "0.0005"},
            {"Alt": "G", "gnomad41_exome_AF_grpmax": "0.008"},
        ]
        selected = [self.rankscore.select_frequency(row, index + 2)[0] for index, row in enumerate(rows)]
        self.assertEqual([0.0005, 0.008], selected)

    def test_nextflow_parameter_and_output_contract(self):
        main_nf = (PIPEVAR_ROOT / "main.nf").read_text()
        config = (PIPEVAR_ROOT / "nextflow.config").read_text()
        readme = (PIPEVAR_ROOT / "README.md").read_text()

        self.assertIn("gnomad_af_ad = null", config)
        self.assertIn("gnomad_af_ar = null", config)
        self.assertIn("gnomad = null", config)
        self.assertIn("0.001", main_nf)
        self.assertIn("0.01", main_nf)
        self.assertIn("cannot be combined with --gnomad_af_ad or --gnomad_af_ar", main_nf)
        self.assertIn("must be less than or equal to --gnomad_af_ar", main_nf)
        self.assertIn("frequency_audit.tsv", readme)

        final_modules = (
            "snp_prio", "multi_snp_prio", "ngs_prio", "multi_ngs_prio",
            "longphase", "multi_longphase", "sv_prio", "multi_sv_prio",
        )
        for module in final_modules:
            module_text = (PIPEVAR_ROOT / "modules" / module / "main.nf").read_text()
            self.assertIn("--gnomad-af-ad", module_text, module)
            self.assertIn("--gnomad-af-ar", module_text, module)
            self.assertIn("--de-novo ${params.denovo_filter}", module_text, module)
            self.assertIn("frequency_audit.tsv", module_text, module)

    def test_vcf_adapters_and_assignment_declare_frequency_metadata(self):
        longphase_root = Path("/home/beoungle/docker_work/longphase")
        for adapter in ("clinvar_vcf_and_txt.sh", "rankscore_vcf_and_txt.sh", "rankvar_vcf_and_tsv.sh"):
            text = (longphase_root / adapter).read_text()
            self.assertIn("ID=PIPEVAR_GNOMAD_AF", text, adapter)
            self.assertIn("ID=PIPEVAR_GNOMAD_SOURCE", text, adapter)
        for assignment in ("assign_dom_or_rec.py", "assign_dom_or_rec_snp_only.py"):
            text = (longphase_root / assignment).read_text()
            self.assertIn('header.info.add("PIPEVAR_GNOMAD_AF"', text, assignment)
            self.assertIn('header.info.add("PIPEVAR_GNOMAD_SOURCE"', text, assignment)


if __name__ == "__main__":
    unittest.main()
