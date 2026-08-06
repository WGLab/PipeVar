import importlib.util
import unittest
from pathlib import Path


MODULE_PATH = Path("/home/beoungle/docker_work/mito_annotation/prioritize_mito_variants.py")
SPEC = importlib.util.spec_from_file_location("prioritize_mito_variants", MODULE_PATH)
MITO = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MITO)


class MitoClinvarScoringTests(unittest.TestCase):
    def test_only_pathogenic_and_likely_pathogenic_receive_support(self):
        accepted = [
            "Pathogenic",
            "Likely_pathogenic",
            "Pathogenic/Likely_pathogenic",
            "Likely-pathogenic,_low-penetrance",
        ]
        rejected = [
            ".",
            "Benign",
            "Likely_benign",
            "Uncertain_significance",
            "Conflicting_classifications_of_pathogenicity",
            "risk_factor",
        ]

        for significance in accepted:
            with self.subTest(significance=significance):
                self.assertEqual(1.0, MITO.score_row({"clinvar_clnsig": significance}))
        for significance in rejected:
            with self.subTest(significance=significance):
                self.assertEqual(0.0, MITO.score_row({"clinvar_clnsig": significance}))


if __name__ == "__main__":
    unittest.main()
