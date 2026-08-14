import csv
import importlib
import sys
import tempfile
import unittest
from pathlib import Path


LONGPHASE_ROOT = Path("/home/beoungle/docker_work/longphase")
PIPEVAR_ROOT = Path(__file__).resolve().parents[1]
if str(LONGPHASE_ROOT) not in sys.path:
    sys.path.insert(0, str(LONGPHASE_ROOT))

prio = importlib.import_module("prio_gene_only")
policy = importlib.import_module("frequency_policy")
phenosv_merge = importlib.import_module("merge_phenosv_event_vcf")


def make_variant(pos, model, af=None, svtype="DEL", gt="0/1", chrom="1", source="phenosv"):
    fields = [
        "Gene=GENE1",
        f"SOURCE={source}",
        f"MODEL={model}",
        "PHENO_SCORE=0.9",
    ]
    if svtype:
        fields.extend([f"SVTYPE={svtype}", f"PIPEVAR_SVTYPE={svtype}"])
    if af is not None:
        fields.extend([
            "COMMON_SV",
            f"COMMON_SV_ID=COMMON_{pos}",
            f"COMMON_SV_AF={af}",
            "COMMON_SV_MATCH=RECIPROCAL_OVERLAP",
            "COMMON_SV_RO=0.8",
        ])
    if not svtype:
        fields.extend(["PIPEVAR_GNOMAD_AF=0.009", "PIPEVAR_GNOMAD_SOURCE=TEST"])
    line = f"{chrom}\t{pos}\tV{pos}\tN\t<{svtype or 'ALT'}>\t.\tPASS\t{';'.join(fields)}\tGT:PS\t{gt}:7"
    return prio.Variant(chrom, str(pos), "N", f"<{svtype or 'ALT'}>", ";".join(fields), "GT:PS", f"{gt}:7", line)


class CommonSvFrequencyPolicyTests(unittest.TestCase):
    def evaluate(self, variant, model, **kwargs):
        return policy.evaluate_frequency(
            variant, model, 0.001, 0.01,
            sv_af_ad=0.005, sv_af_ar=0.01, **kwargs,
        )

    def test_exact_boundaries_pass(self):
        self.assertEqual("PASS", self.evaluate(make_variant(10, "Dominant", "0.005"), "Dominant")[0])
        self.assertEqual("PASS", self.evaluate(make_variant(11, "Recessive", "0.01"), "Recessive")[0])

    def test_intermediate_frequency_prunes_only_dominant_lane(self):
        variant = make_variant(12, "Dominant;Recessive", "0.008", gt="1/1")
        scenarios = prio.get_all_gene_scenarios(
            [variant], common_sv_af_ad=0.005, common_sv_af_ar=0.01,
        )
        self.assertEqual(["PhenoSV_AR_Hom"], [item["cat"] for item in scenarios])

    def test_xlinked_and_de_novo_threshold_selection(self):
        xld = make_variant(121, "XLD", "0.008", chrom="X")
        xlr = make_variant(122, "XLR", "0.008", chrom="X", gt="1/1")
        de_novo = make_variant(123, "Recessive", "0.008", gt="1/1")
        self.assertEqual("FAIL", self.evaluate(xld, "XLD")[0])
        self.assertEqual("PASS", self.evaluate(xlr, "XLR")[0])
        self.assertEqual("FAIL", self.evaluate(de_novo, "Recessive", de_novo=True)[0])

    def test_unknown_duplication_uses_dominant_ceiling(self):
        variant = make_variant(13, "Duplication", "0.008", svtype="DUP")
        self.assertEqual([], prio.get_all_gene_scenarios([variant]))

    def test_clinvar_only_unknown_moi_uses_dominant_ceiling(self):
        passing = make_variant(131, "ClinVar_SV", "0.005", source="clinvar")
        failing = make_variant(132, "ClinVar_SV", "0.008", source="clinvar")
        self.assertEqual(
            ["ClinVar_AD"],
            [item["cat"] for item in prio.get_all_gene_scenarios([passing])],
        )
        self.assertEqual([], prio.get_all_gene_scenarios([failing]))

    def test_disabled_filter_and_bnd_are_not_evaluated(self):
        disabled = self.evaluate(
            make_variant(14, "Dominant", "0.8"), "Dominant",
            common_sv_evaluated=False,
        )
        bnd = self.evaluate(make_variant(15, "Dominant", "0.8", svtype="BND"), "Dominant")
        self.assertEqual(("PASS", None, "COMMON_SV_NOT_EVALUATED"), disabled)
        self.assertEqual(("PASS", None, "SVTYPE_NOT_EVALUATED"), bnd)

    def test_mitochondrial_sv_is_exempt(self):
        self.assertEqual(
            ("PASS", None, "MITO_EXEMPT"),
            self.evaluate(make_variant(16, "Dominant", "0.8", chrom="MT"), "Dominant"),
        )

    def test_mixed_snv_sv_comphet_members_pass_independently(self):
        snv = make_variant(20, "Recessive", svtype="", gt="1|0")
        sv = make_variant(21, "Recessive", "0.008", gt="0|1")
        scenarios = prio.get_all_gene_scenarios([snv, sv])
        self.assertEqual(1, len(scenarios))
        self.assertIn("_AR+", scenarios[0]["cat"])

    def test_frequency_failure_remains_in_audit(self):
        variant = make_variant(30, "Dominant", "0.008")
        with tempfile.TemporaryDirectory() as tmpdir:
            audit = Path(tmpdir) / "audit.tsv"
            prio.write_frequency_audit(
                [variant], audit, 0.001, 0.01,
                common_sv_af_ad=0.005, common_sv_af_ar=0.01,
            )
            with audit.open(newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual("FAIL", rows[0]["decision"])
        self.assertEqual("0.008", rows[0]["selected_af"])
        self.assertEqual("PIPEVAR_COMMON_SV_COHORT", rows[0]["af_source"])

    def test_nextflow_parameter_and_container_contract(self):
        config = (PIPEVAR_ROOT / "nextflow.config").read_text()
        main = (PIPEVAR_ROOT / "main.nf").read_text()
        self.assertIn('pipeline_version = "0.5.0"', config)
        self.assertIn("common_sv_af = null", config)
        self.assertIn("common_sv_af_ad = null", config)
        self.assertIn("common_sv_af_ar = null", config)
        self.assertIn("cannot be combined with --common_sv_af_ad or --common_sv_af_ar", main)
        self.assertIn("must be less than or equal to --common_sv_af_ar", main)

        for module in ("common_sv_filter", "multi_common_sv_filter"):
            text = (PIPEVAR_ROOT / "modules" / module / "main.nf").read_text()
            self.assertIn("common_sv_filter_0.5", text)
            self.assertIn("--common-af-ad", text)
            self.assertIn("--common-af-ar", text)

        for module in (
            "snp_prio", "multi_snp_prio", "ngs_prio", "multi_ngs_prio",
            "longphase", "multi_longphase", "sv_prio", "multi_sv_prio",
        ):
            text = (PIPEVAR_ROOT / "modules" / module / "main.nf").read_text()
            self.assertIn("longphase_0.4.0", text)
            self.assertIn("--common-sv-af-ad ${params.common_sv_af_ad}", text)
            self.assertIn("--common-sv-af-ar ${params.common_sv_af_ar}", text)
            self.assertIn("--common-sv-filter ${params.common_sv_filter}", text)

    def test_provenance_fields_are_carried_by_merge_and_assignment(self):
        merge_source = (LONGPHASE_ROOT / "merge_phenosv_event_vcf.py").read_text()
        for field in (
            "COMMON_SV_ID", "COMMON_SV_AF", "COMMON_SV_MATCH",
            "COMMON_SV_RO", "COMMON_SV_INS_IDENTITY",
        ):
            self.assertIn(field, merge_source)
            for assignment in ("assign_dom_or_rec.py", "assign_dom_or_rec_sv_only.py"):
                self.assertIn(field, (LONGPHASE_ROOT / assignment).read_text())

    def test_phenosv_merge_copies_common_sv_provenance_from_annotation(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            events = root / "events.tsv"
            base = root / "base.vcf"
            annotation = root / "annotation.vcf"
            output = root / "output.vcf"
            events.write_text(
                "SV_ID\tPHENOSV_EVENT_ID\tCHROM\tSTART\tEND\tSVTYPE\tPATHOGENICITY\tPHEN2GENE\tPHENOSV_SCORE\tPHENOSV_TYPE\tPHENOSV_GENE\tPHENOSV_GENE_SCORE\n"
                "SV1\tEVENT1\t1\t100\t200\tDEL\t0.9\t0.8\t0.85\texonic\tGENE1\t0.7\n"
            )
            header = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            base.write_text(header + "1\t100\tSV1\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=200\n")
            annotation.write_text(
                header
                + "1\t100\tSV1\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=200;Gene.refGene=GENE1;COMMON_SV;COMMON_SV_ID=DB1;COMMON_SV_AF=0.008;COMMON_SV_MATCH=RECIPROCAL_OVERLAP;COMMON_SV_RO=0.8\n"
            )
            phenosv_merge.merge(events, base, annotation, output)
            merged = output.read_text()
        self.assertIn("COMMON_SV_ID=DB1", merged)
        self.assertIn("COMMON_SV_AF=0.008", merged)
        self.assertIn("COMMON_SV_MATCH=RECIPROCAL_OVERLAP", merged)
        self.assertIn("COMMON_SV_RO=0.8", merged)


if __name__ == "__main__":
    unittest.main()
