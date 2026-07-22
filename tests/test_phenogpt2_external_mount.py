import unittest
from pathlib import Path


PIPEVAR = Path(__file__).resolve().parents[1]
DOCKER_CONTEXT = Path("/home/beoungle/docker_work/phenogpt2")


class PhenoGPT2ExternalMountTests(unittest.TestCase):
    def test_docker_context_excludes_models_and_generated_caches(self):
        dockerfile = (DOCKER_CONTEXT / "Dockerfile").read_text()
        dockerignore = (DOCKER_CONTEXT / ".dockerignore").read_text().splitlines()
        self.assertNotIn("COPY new_model", dockerfile)
        self.assertIn("new_model/", dockerignore)
        self.assertIn(".cache/", dockerignore)
        self.assertIn("models--*/", dockerignore)

    def test_modules_use_same_fixed_image_and_internal_paths(self):
        modules = [
            (PIPEVAR / "modules/phenogpt2/main.nf").read_text(),
            (PIPEVAR / "modules/multi_phenogpt2/main.nf").read_text(),
        ]
        for module in modules:
            self.assertIn("beoungl/docker_test:phenogpt2_0.2", module)
            self.assertIn("/opt/phenogpt2/models/phenogpt2", module)
            self.assertIn("PHENOGPT2_MODEL_FINGERPRINT", module)
        self.assertNotIn("params.phenogpt2_container", (PIPEVAR / "nextflow.config").read_text())

    def test_mounts_are_phenogpt2_process_specific(self):
        config = (PIPEVAR / "nextflow.config").read_text()
        deepvariant_block = config.split("withName: 'deepvariant|multi_deepvariant'", 1)[1].split("}", 1)[0]
        self.assertNotIn("phenogpt2_model_host_path", deepvariant_block)
        self.assertIn(":/opt/phenogpt2/models/phenogpt2:ro", config)
        self.assertIn(":/opt/phenogpt2/cache:rw", config)

    def test_legacy_csv_defaults_to_clinical_notes(self):
        main = (PIPEVAR / "main.nf").read_text()
        self.assertIn(
            "!(params.note != null && clean_note == 'no')",
            main,
        )


if __name__ == "__main__":
    unittest.main()
