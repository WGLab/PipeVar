import csv
import subprocess
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path("/home/beoungle/docker_work/rankscore/clinvar.sh")


class RankscoreClinvarFilterTests(unittest.TestCase):
    def test_clnsig_filter_is_header_indexed_and_preserves_raw_rows(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            annovar = root / "sample.multianno.txt"
            phen2gene = root / "phen2gene.tsv"
            prefix = root / "sample"
            header = [
                "Marker", "Gene.refGene", "CLNSIG", "Predictor", "Func.refGene",
                "gnomad41_exome_AF_raw", "gnomad41_genome_AF_raw", "Otherinfo12", "Otherinfo13",
            ]
            assertions = [
                ("p", "Pathogenic", "."),
                ("lp", "likely-pathogenic", "."),
                ("compound", "Pathogenic/Likely_pathogenic", "."),
                ("low_penetrance", "Likely_pathogenic,_low_penetrance", "."),
                ("conflicting", "Conflicting_classifications_of_pathogenicity", "."),
                ("vus", "Uncertain_significance", "Pathogenic"),
                ("benign", "Benign", "."),
                ("risk", "risk_factor", "."),
            ]
            with annovar.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                writer.writerow(header)
                for marker, clnsig, predictor in assertions:
                    writer.writerow([marker, "GENE1", clnsig, predictor, "intronic", ".", ".", "GT", "0/1"])
            phen2gene.write_text("Rank\tGene\n1\tGENE1\n", encoding="utf-8")

            subprocess.run(
                [
                    "bash", str(SCRIPT), str(annovar), str(phen2gene), str(prefix),
                    "0.0001", "0.5", "1", "20", "15", "REVEL",
                ],
                text=True,
                capture_output=True,
            )

            rows = list(csv.DictReader(
                Path(f"{prefix}.clinvar.txt").read_text(encoding="utf-8").splitlines(),
                delimiter="\t",
            ))
            self.assertEqual(
                {"p", "lp", "compound", "low_penetrance"},
                {row["Marker"] for row in rows},
            )
            self.assertEqual("likely-pathogenic", next(row for row in rows if row["Marker"] == "lp")["CLNSIG"])


if __name__ == "__main__":
    unittest.main()
