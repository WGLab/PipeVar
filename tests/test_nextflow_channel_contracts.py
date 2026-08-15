from __future__ import annotations

import shutil
import subprocess
import tempfile
import textwrap
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
NEXTFLOW = shutil.which("nextflow")


class FailFastFanInContractTests(unittest.TestCase):
    def test_required_batch_fan_ins_are_fail_fast(self):
        required_joins = {
            "subworkflows/input_csv_annotated_snv_sv/main.nf": (
                "validated_annovar_txt.join(phen2gene_result, failOnMismatch: true, failOnDuplicate: true)",
                "join_annovar_phen2gene.join(hpo_paths, failOnMismatch: true, failOnDuplicate: true)",
                "survivor_result.join(hpo_paths, failOnMismatch: true, failOnDuplicate: true)",
                "phenosv_result.join(validated_annovar_vcf, failOnMismatch: true, failOnDuplicate: true)",
                "phenosv_annovar_snv.join(annovar_sv_for_downstream, failOnMismatch: true, failOnDuplicate: true)",
                "sv_join.join(rankscore_result, failOnMismatch: true, failOnDuplicate: true)",
                "rankscore_join.join(rankvar_result, failOnMismatch: true, failOnDuplicate: true)",
                "rankvar_join.join(hpo_with_age, failOnMismatch: true, failOnDuplicate: true)",
            ),
            "subworkflows/input_csv_annotated_all_ngs/main.nf": (
                "validated_annovar_txt.join(phen2gene_result, failOnMismatch: true, failOnDuplicate: true)",
                "join_annovar_phen2gene.join(hpo_paths, failOnMismatch: true, failOnDuplicate: true)",
                "survivor_result.join(hpo_paths, failOnMismatch: true, failOnDuplicate: true)",
                "phenosv_result.join(validated_annovar_vcf, failOnMismatch: true, failOnDuplicate: true)",
                "phenosv_annovar_snv.join(annovar_sv_for_downstream, failOnMismatch: true, failOnDuplicate: true)",
                "sv_join.join(rankscore_result, failOnMismatch: true, failOnDuplicate: true)",
                "rankscore_join.join(rankvar_result, failOnMismatch: true, failOnDuplicate: true)",
                "rankvar_join.join(hpo_with_age, failOnMismatch: true, failOnDuplicate: true)",
            ),
            "subworkflows/input_csv_annotated_snv_called_sv_ngs/main.nf": (
                "validated_annovar_txt.join(phen2gene_result, failOnMismatch: true, failOnDuplicate: true)",
                "join_annovar_phen2gene.join(hpo_paths, failOnMismatch: true, failOnDuplicate: true)",
                "survivor_result.join(hpo_paths, failOnMismatch: true, failOnDuplicate: true)",
                "phenosv_result.join(validated_annovar_vcf, failOnMismatch: true, failOnDuplicate: true)",
                "phenosv_annovar_snv.join(annovar_sv_for_downstream, failOnMismatch: true, failOnDuplicate: true)",
                "sv_join.join(rankscore_result, failOnMismatch: true, failOnDuplicate: true)",
                "rankscore_join.join(rankvar_result, failOnMismatch: true, failOnDuplicate: true)",
                "rankvar_join.join(hpo_with_age, failOnMismatch: true, failOnDuplicate: true)",
            ),
            "subworkflows/input_csv_alignment_ngs_sv/main.nf": (
                "survivor_result.join(input_bam_no_bam, failOnMismatch: true, failOnDuplicate: true)",
                "phenosv_result.join(annovar_sv_for_downstream, failOnMismatch: true, failOnDuplicate: true)",
                "sv_prio_input.join(input_bam_hpo_age, failOnMismatch: true, failOnDuplicate: true)",
            ),
            "subworkflows/input_csv_alignment_long_sv/main.nf": (
                "survivor_result.join(input_bam_no_bam, failOnMismatch: true, failOnDuplicate: true)",
                "phenosv_result.join(annovar_sv_for_downstream, failOnMismatch: true, failOnDuplicate: true)",
                "sv_prio_input.join(input_bam_hpo_age, failOnMismatch: true, failOnDuplicate: true)",
            ),
            "subworkflows/input_csv_alignment_vcf_sv/main.nf": (
                "survivor_result.join(input_vcf_no_vcf, failOnMismatch: true, failOnDuplicate: true)",
                "phenosv_result.join(annovar_sv_for_downstream, failOnMismatch: true, failOnDuplicate: true)",
                "sv_prio_input.join(input_vcf_hpo_age, failOnMismatch: true, failOnDuplicate: true)",
            ),
            "subworkflows/input_csv_alignment_all_ngs/main.nf": (
                "survivor_result.join(input_bam_no_bam, failOnMismatch: true, failOnDuplicate: true)",
                "phenosv_result.join(annovar_result_vcf, failOnMismatch: true, failOnDuplicate: true)",
                "phenosv_annovar_snv.join(annovar_sv_for_downstream, failOnMismatch: true, failOnDuplicate: true)",
                "sv_join.join(rankscore_result, failOnMismatch: true, failOnDuplicate: true)",
                "rankscore_join.join(rankvar_result, failOnMismatch: true, failOnDuplicate: true)",
                "rankvar_join.join(input_bam_hpo_age, failOnMismatch: true, failOnDuplicate: true)",
            ),
            "subworkflows/input_csv_alignment_all_ngs_light/main.nf": (
                "survivor_result.join(input_bam_no_bam, failOnMismatch: true, failOnDuplicate: true)",
                "phenosv_result.join(annovar_result_vcf, failOnMismatch: true, failOnDuplicate: true)",
                "phenosv_annovar_snv.join(annovar_sv_result, failOnMismatch: true, failOnDuplicate: true)",
                "sv_join.join(rankscore_result, failOnMismatch: true, failOnDuplicate: true)",
                "rankscore_join.join(rankvar_result, failOnMismatch: true, failOnDuplicate: true)",
                "rankvar_join.join(input_bam_hpo_age, failOnMismatch: true, failOnDuplicate: true)",
            ),
            "subworkflows/input_csv_alignment_all_longphase/main.nf": (
                "survivor_result.join(input_bam_no_bam, failOnMismatch: true, failOnDuplicate: true)",
                "rankscore_result.join(join_vcf_bam_phenosv, failOnMismatch: true, failOnDuplicate: true)",
                "rankvar_result.join(join_vcf_bam_rankscore, failOnMismatch: true, failOnDuplicate: true)",
                "join_vcf_bam_rankvar.join(input_bam_hpo_age, failOnMismatch: true, failOnDuplicate: true)",
            ),
            "subworkflows/input_csv_alignment_all_light_longphase/main.nf": (
                "survivor_result.join(input_bam_no_bam, failOnMismatch: true, failOnDuplicate: true)",
                "rankscore_result.join(join_vcf_bam_phenosv, failOnMismatch: true, failOnDuplicate: true)",
                "rankvar_result.join(join_vcf_bam_rankscore, failOnMismatch: true, failOnDuplicate: true)",
                "join_vcf_bam_rankvar.join(input_bam_hpo_age, failOnMismatch: true, failOnDuplicate: true)",
            ),
        }

        for relative_path, joins in required_joins.items():
            source = (ROOT / relative_path).read_text(encoding="utf-8")
            for join in joins:
                with self.subTest(path=relative_path, join=join):
                    self.assertIn(join, source)


@unittest.skipUnless(NEXTFLOW, "Nextflow is not installed")
class NextflowRuntimeContractTests(unittest.TestCase):
    def run_nextflow(self, script: str, temp_dir: Path, *args: str) -> subprocess.CompletedProcess[str]:
        config = temp_dir / "empty.config"
        config.write_text("nextflow.enable.dsl=2\n", encoding="utf-8")
        workflow = temp_dir / "main.nf"
        workflow.write_text(textwrap.dedent(script), encoding="utf-8")
        return subprocess.run(
            [
                NEXTFLOW,
                "-C",
                str(config),
                "run",
                str(workflow),
                "-ansi-log",
                "false",
                "-work-dir",
                str(temp_dir / "work"),
                *args,
            ],
            cwd=ROOT,
            text=True,
            capture_output=True,
            check=False,
        )

    def test_validator_emits_two_tuples_to_downstream_process(self):
        with tempfile.TemporaryDirectory(dir=ROOT) as temp_name:
            temp_dir = Path(temp_name)
            results = temp_dir / "results"
            results.mkdir()
            inputs = []
            for sample in ("case1", "case2"):
                txt = temp_dir / f"{sample}.hg38_multianno.txt"
                vcf = temp_dir / f"{sample}.hg38_multianno.vcf"
                txt.write_text("Chr\tStart\n", encoding="utf-8")
                vcf.write_text("##fileformat=VCFv4.2\n", encoding="utf-8")
                inputs.extend((f"--{sample}_txt", str(txt), f"--{sample}_vcf", str(vcf)))

            script = '''
                nextflow.enable.dsl=2

                include { validate_preannotated_annovar_pair } from '../modules/validate_preannotated_annovar_pair'

                process validation_sentinel {
                    publishDir params.results_dir, mode: 'copy'

                    input:
                    tuple val(out_prefix), path(annovar_txt), path(annovar_vcf)

                    output:
                    path "${out_prefix}.seen"

                    stub:
                    """
                    printf '%s\\n' '$out_prefix' > ${out_prefix}.seen
                    """
                }

                workflow {
                    input_pairs = channel.of(
                        tuple('case1', file(params.case1_txt), file(params.case1_vcf)),
                        tuple('case2', file(params.case2_txt), file(params.case2_vcf))
                    )
                    validated_pairs = validate_preannotated_annovar_pair(input_pairs)
                    validation_sentinel(validated_pairs)
                }
            '''
            result = self.run_nextflow(
                script,
                temp_dir,
                "-stub-run",
                "--results_dir",
                str(results),
                *inputs,
            )

            self.assertEqual(0, result.returncode, result.stdout + result.stderr)
            self.assertEqual(
                {"case1.seen", "case2.seen"},
                {path.name for path in results.glob("*.seen")},
            )

    def test_join_rejects_missing_and_duplicate_keys(self):
        cases = {
            "missing": "right = channel.of(tuple('case1', 'R1'))",
            "duplicate": "right = channel.of(tuple('case1', 'R1'), tuple('case1', 'R2'), tuple('case2', 'R3'))",
        }
        for name, right_channel in cases.items():
            with self.subTest(name=name), tempfile.TemporaryDirectory(dir=ROOT) as temp_name:
                script = f"""
                    nextflow.enable.dsl=2

                    workflow {{
                        left = channel.of(tuple('case1', 'L1'), tuple('case2', 'L2'))
                        {right_channel}
                        left.join(right, failOnMismatch: true, failOnDuplicate: true).view()
                    }}
                """
                result = self.run_nextflow(script, Path(temp_name))
                self.assertNotEqual(0, result.returncode, result.stdout + result.stderr)


if __name__ == "__main__":
    unittest.main()
