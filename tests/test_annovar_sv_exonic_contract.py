import re
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
SINGLE_MODULE = REPO / "modules/annovar_sv/main.nf"
BATCH_MODULE = REPO / "modules/multi_annovar_sv/main.nf"
SV_CALLSITE_FILES = tuple(
    sorted(
        path
        for path in (REPO / "subworkflows").glob("*/main.nf")
        if re.search(r"(?:^|\W)(?:ANNOVAR_SV|multi_annovar_sv)\(", path.read_text(encoding="utf-8"))
    )
)


class AnnovarSvExonicContractTests(unittest.TestCase):
    def test_modules_filter_raw_multianno_vcf_without_phen2gene_or_bed(self):
        single = SINGLE_MODULE.read_text(encoding="utf-8")
        batch = BATCH_MODULE.read_text(encoding="utf-8")

        self.assertIn("path vcf", single)
        self.assertIn("val out_prefix", single)
        self.assertIn("val sv_annotation_mode", single)
        self.assertIn("tuple val(out_prefix), path(vcf), val(sv_annotation_mode)", batch)
        for source in (single, batch):
            self.assertIn("bcftools view -i 'INFO/Func.refGene=\"exonic\"'", source)
            self.assertNotIn("phen2gene", source.lower())
            self.assertNotIn("bed_file", source)
            self.assertNotIn("-bedfile", source)

    def test_all_sv_call_sites_use_the_three_field_contract(self):
        self.assertEqual(19, len(SV_CALLSITE_FILES))
        sources = "\n".join(path.read_text(encoding="utf-8") for path in SV_CALLSITE_FILES)

        self.assertNotRegex(sources, r"ANNOVAR_SV\([^\n]*(?:Phen2gene|phen2gene|bed_file|phen2_gene_bed)")
        self.assertNotRegex(
            sources,
            r"tuple\([^\n]*vcf_file[^\n]*(?:phen2gene_file|bed_file)[^\n]*(?:called|preannotated)",
        )
        self.assertEqual(9, len(re.findall(r"(?<!multi_)ANNOVAR_SV\([^\n]+\)", sources)))
        self.assertEqual(10, len(re.findall(r"multi_annovar_sv\([^\n]+\)", sources)))

    def test_standalone_sv_workflows_do_not_run_phen2gene(self):
        standalone = (
            "subworkflows/single_alignment_ngs_sv/main.nf",
            "subworkflows/single_alignment_long_sv/main.nf",
            "subworkflows/single_alignment_vcf_sv/main.nf",
            "subworkflows/input_csv_alignment_ngs_sv/main.nf",
            "subworkflows/input_csv_alignment_long_sv/main.nf",
            "subworkflows/input_csv_alignment_vcf_sv/main.nf",
        )
        for relative in standalone:
            source = (REPO / relative).read_text(encoding="utf-8")
            self.assertNotRegex(source, r"include \{ (?:multi_)?[Pp]hen2gene")
            self.assertNotRegex(source, r"(?:Phen2gene|multi_phen2gene)\(")

    @unittest.skipUnless(shutil.which("bcftools"), "bcftools is provided by the truvari runtime image")
    def test_exact_filter_retains_multigene_exonic_records_in_order(self):
        fixture = """\
##fileformat=VCFv4.2
##contig=<ID=1>
##INFO=<ID=Func.refGene,Number=.,Type=String,Description="ANNOVAR functional category">
##INFO=<ID=Gene.refGene,Number=.,Type=String,Description="ANNOVAR gene annotation">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="SV type">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
1\t100\tmulti_gene_del\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;Func.refGene=exonic;Gene.refGene=GENEA,GENEB,GENEC
1\t300\tunmatched_gene_del\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;Func.refGene=exonic;Gene.refGene=OTHER
1\t500\tintronic_del\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;Func.refGene=intronic;Gene.refGene=GENED
1\t700\tsplicing_del\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;Func.refGene=splicing;Gene.refGene=GENEE
"""
        with tempfile.TemporaryDirectory() as tmpdir:
            input_vcf = Path(tmpdir) / "input.vcf"
            output_vcf = Path(tmpdir) / "output.vcf"
            input_vcf.write_text(fixture, encoding="utf-8")
            subprocess.run(
                [
                    "bcftools",
                    "view",
                    "-i",
                    'INFO/Func.refGene="exonic"',
                    "-Ov",
                    "-o",
                    str(output_vcf),
                    str(input_vcf),
                ],
                check=True,
            )

            records = [line for line in output_vcf.read_text(encoding="utf-8").splitlines() if not line.startswith("#")]
            self.assertEqual(["multi_gene_del", "unmatched_gene_del"], [line.split("\t")[2] for line in records])
            self.assertIn("Gene.refGene=GENEA,GENEB,GENEC", records[0])


if __name__ == "__main__":
    unittest.main()
