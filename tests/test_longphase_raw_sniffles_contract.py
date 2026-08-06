import re
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
DOCKER_WORK = REPO.parent / "docker_work"


def read(relative_path):
    return (REPO / relative_path).read_text(encoding="utf-8")


class RawSnifflesLongPhaseContractTests(unittest.TestCase):
    def test_sniffles_modules_publish_rnames_with_versioned_image(self):
        for module in ("modules/sniffles/main.nf", "modules/multi_sniffles/main.nf"):
            source = read(module)
            self.assertIn("beoungl/docker_test:sniffles_0.1", source)
            self.assertIn("--output-rnames", source)
            self.assertIn(".sniffles.vcf", source)

        dockerfile = (DOCKER_WORK / "sniffles/Dockerfile").read_text(encoding="utf-8")
        self.assertIn("sniffles_0.1", dockerfile)
        self.assertIn("sniffles:2.8.0--", dockerfile)

    def test_prioritization_and_phasing_modules_adopt_the_new_image(self):
        expected = {
            REPO / "modules/longphase/main.nf",
            REPO / "modules/multi_longphase/main.nf",
            REPO / "modules/ngs_prio/main.nf",
            REPO / "modules/multi_ngs_prio/main.nf",
            REPO / "modules/snp_prio/main.nf",
            REPO / "modules/multi_snp_prio/main.nf",
            REPO / "modules/sv_prio/main.nf",
            REPO / "modules/multi_sv_prio/main.nf",
        }
        actual = {
            path
            for path in (REPO / "modules").glob("*/main.nf")
            if "longphase_0.2.32" in path.read_text(encoding="utf-8")
        }
        self.assertEqual(expected, actual)
        for path in (REPO / "modules/longphase/main.nf", REPO / "modules/multi_longphase/main.nf"):
            source = path.read_text(encoding="utf-8")
            self.assertRegex(source, r"path\(sv_phase_vcf\)|path sv_phase_vcf")
            self.assertIn("--sv-file=$sv_phase_vcf", source)
            self.assertIn("task.ext.args != null ? task.ext.args : '--ont'", source)

        config = read("nextflow.config")
        self.assertRegex(config, r"clean_type == 'ont'[\s\S]+?ext\.args='--ont'")
        self.assertRegex(config, r"clean_type == 'pacbio'[\s\S]+?ext\.args=''")

    def test_single_sample_phases_raw_sniffles_but_curates_evidence(self):
        source = read("subworkflows/single_alignment_all_longphase/main.nf")
        self.assertIn("ANNOVAR_SV(sniffles.out", source)
        self.assertIn("SURVIVOR(annovar_sv_for_downstream", source)
        self.assertIn("PhenoSV(SURVIVOR.out", source)
        self.assertRegex(
            source,
            r"longphase\(bam,ANNOVAR\.out\.vcf_output,sniffles\.out,annovar_sv_for_downstream,PhenoSV\.out,",
        )

    def test_batch_phasing_is_proband_keyed_and_annotation_remains_curated(self):
        source = read("subworkflows/input_csv_alignment_all_longphase/main.nf")
        self.assertIn("sv_for_annotation=denovo_sv_result.records", source)
        self.assertIn(
            "sniffles_for_phasing=sniffles_result.join(proband_keys, failOnDuplicate: true)",
            source,
        )
        self.assertIn("annovar_sv_result=multi_annovar_sv(sniffles_result_annovar)", source)
        self.assertIn("survivor_result=multi_survivor(annovar_sv_for_downstream)", source)
        self.assertIn("phenosv_result=multi_phenosv(phenosv_input)", source)
        self.assertIn(
            "sniffles_for_phasing.join(annovar_sv_for_downstream, failOnMismatch: true, failOnDuplicate: true).join(join_vcf_bam, failOnMismatch: true, failOnDuplicate: true)",
            source,
        )
        self.assertIn(
            "phenosv_result.join(join_vcf_bam_sv, failOnMismatch: true, failOnDuplicate: true)",
            source,
        )

    def test_light_longphase_routes_also_phase_raw_sniffles(self):
        single = read("subworkflows/single_alignment_all_light_longphase/main.nf")
        batch = read("subworkflows/input_csv_alignment_all_light_longphase/main.nf")
        self.assertRegex(single, r"longphase\(bam,ANNOVAR\.out\.vcf_output,sniffles\.out,ANNOVAR_SV\.out,PhenoSV\.out,")
        self.assertIn("join_vcf_bam_sv=sniffles_result.join(", batch)

    def test_cnvpytor_and_longread_merge_are_inactive_but_modules_remain(self):
        active_paths = [REPO / "main.nf", REPO / "nextflow.config", REPO / "README.md"]
        active_paths.extend((REPO / "docs").rglob("*"))
        active_paths.extend((REPO / "subworkflows").rglob("*.nf"))
        for path in active_paths:
            if not path.is_file():
                continue
            source = path.read_text(encoding="utf-8")
            self.assertIsNone(re.search(r"cnvpytor|merge_longread_sv_callers", source, re.I), str(path))

        retained = (
            "modules/cnvpytor/main.nf",
            "modules/multi_cnvpytor/main.nf",
            "modules/merge_longread_sv_callers/main.nf",
            "modules/multi_merge_longread_sv_callers/main.nf",
        )
        for relative_path in retained:
            self.assertTrue((REPO / relative_path).is_file(), relative_path)
        self.assertIn("quay.io/biocontainers/cnvpytor:1.3.2--pyhdfd78af_0", read(retained[0]))
        self.assertIn("quay.io/biocontainers/cnvpytor:1.3.2--pyhdfd78af_0", read(retained[1]))

    def test_longphase_container_preserves_ps_and_requires_same_phase_set(self):
        assign = (DOCKER_WORK / "longphase/assign_dom_or_rec.py").read_text(encoding="utf-8")
        assign_sv = (DOCKER_WORK / "longphase/assign_dom_or_rec_sv_only.py").read_text(encoding="utf-8")
        prioritizer = (DOCKER_WORK / "longphase/prio_gene_only.py").read_text(encoding="utf-8")
        assign_snp = (DOCKER_WORK / "longphase/assign_dom_or_rec_snp_only.py").read_text(encoding="utf-8")
        for source in (assign, assign_sv, assign_snp):
            self.assertIn('header.formats.add("PS"', source)
            self.assertRegex(source, r"row(?:\[|\.get\()['\"]ps['\"]")
            self.assertRegex(source, r"'id': (?:parts\[2\]|record_id)")
            self.assertIn("rec.id = str(row['id'])", source)
        self.assertIn("ps1 is not None and ps1 == ps2", prioritizer)
        self.assertIn("Every other non-trans heterozygous pairing is unresolved", prioritizer)


if __name__ == "__main__":
    unittest.main()
