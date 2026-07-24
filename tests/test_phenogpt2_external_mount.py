import unittest
from pathlib import Path


PIPEVAR = Path(__file__).resolve().parents[1]
DOCKER_CONTEXT = Path("/home/beoungle/docker_work/phenogpt2")


class PhenoGPT2ExternalMountTests(unittest.TestCase):
    def test_docker_context_excludes_models_and_generated_caches(self):
        dockerfile = (DOCKER_CONTEXT / "Dockerfile").read_text()
        dockerignore = (DOCKER_CONTEXT / ".dockerignore").read_text().splitlines()
        self.assertNotIn("COPY new_model", dockerfile)
        self.assertNotIn("phenogpt2-runtime.patch", dockerfile)
        self.assertIn("new_model/", dockerignore)
        self.assertIn(".cache/", dockerignore)
        self.assertIn("models--*/", dockerignore)

    def test_dockerfile_has_pinned_slim_targets(self):
        dockerfile = (DOCKER_CONTEXT / "Dockerfile").read_text()
        self.assertIn("sha256:09d8951b943dee03cf8fc841b6ea1f201ad33f82f76567171394853c0f494054", dockerfile)
        self.assertIn("sha256:a5dc56795348ea8fe46305612afaf7f9882cb438b08a582d4e5bd762369d7922", dockerfile)
        self.assertIn("AS artifact-slim", dockerfile)
        self.assertIn("AS dependency-slim", dockerfile)
        self.assertIn("AS runtime-slim", dockerfile)
        self.assertIn("MINICONDA_SHA256=a098a5b1581d8fd078c430b82e27106602223e335efef708a124e723814d120c", dockerfile)
        self.assertIn("PHENOGPT2_REF=68922bc131abecc3e3becd0e5ae927d015e4b0e4", dockerfile)

    def test_runtime_requirements_keep_negation_and_drop_training_packages(self):
        requirements = (DOCKER_CONTEXT / "requirements-runtime.txt").read_text().lower()
        lock = (DOCKER_CONTEXT / "requirements-runtime.lock").read_text().lower()
        dockerfile = (DOCKER_CONTEXT / "Dockerfile").read_text()
        for package in (
            "sentence-transformers",
            "rapidfuzz",
            "scipy",
            "scikit-learn",
            "joblib",
            "spacy",
            "vllm",
            "xformers",
        ):
            self.assertIn(package, requirements)
        for package in (
            "matplotlib",
            "seaborn",
            "peft",
            "trl==",
            "deepspeed",
            "datasets",
            "bitsandbytes",
            "pyarrow",
            "nltk",
            "networkx",
            "obonet",
        ):
            self.assertNotIn(package, requirements)
            if package != "networkx":
                self.assertNotIn(package, lock)
        self.assertIn("-r /tmp/requirements-runtime.lock", dockerfile)

    def test_modules_use_same_fixed_image_and_internal_paths(self):
        modules = [
            (PIPEVAR / "modules/phenogpt2/main.nf").read_text(),
            (PIPEVAR / "modules/multi_phenogpt2/main.nf").read_text(),
        ]
        for module in modules:
            self.assertIn("beoungl/docker_test:phenogpt2_0.2", module)
            self.assertIn("/opt/phenogpt2/models/phenogpt2", module)
            self.assertIn("/opt/phenogpt2/models/negation", module)
            self.assertIn("/opt/phenogpt2/models/embedding", module)
            self.assertNotIn("PHENOGPT2_MODEL_FINGERPRINT", module)
            self.assertNotIn("PHENOGPT2_COMBINED_MODEL_FINGERPRINT", module)
        self.assertNotIn("params.phenogpt2_container", (PIPEVAR / "nextflow.config").read_text())

    def test_mounts_are_phenogpt2_process_specific(self):
        config = (PIPEVAR / "nextflow.config").read_text()
        deepvariant_block = config.split("withName: 'deepvariant|multi_deepvariant'", 1)[1].split("}", 1)[0]
        self.assertNotIn("phenogpt2_model_host_path", deepvariant_block)
        self.assertIn(":/opt/phenogpt2/models/phenogpt2:ro", config)
        self.assertIn(":/opt/phenogpt2/models/negation:ro", config)
        self.assertIn(":/opt/phenogpt2/models/embedding:ro", config)
        self.assertIn(":/opt/phenogpt2/cache:rw", config)

    def test_negation_is_enabled_without_upstream_embedding_argument(self):
        wrapper = (DOCKER_CONTEXT / "run_phenogpt2_to_hpo").read_text()
        helper = (DOCKER_CONTEXT / "phenogpt2_to_hpo.py").read_text()
        main = (PIPEVAR / "main.nf").read_text()
        self.assertIn("-negation_model \"$negation_model_dir\"", wrapper)
        self.assertNotIn("-embedding_model", wrapper)
        self.assertIn("prepare-embedding-cache", wrapper)
        self.assertIn("models--Qwen--Qwen3-Embedding-0.6B", helper)
        self.assertIn("PhenoGPT2 requires exactly one visible CUDA GPU", wrapper)
        self.assertNotIn("negation yes is not supported", main)
        self.assertIn("phenogpt2_negation_model_host_path", main)
        self.assertIn("phenogpt2_embedding_model_host_path", main)
        self.assertNotIn("SHA256SUMS", main)
        self.assertNotIn("MessageDigest", main)
        self.assertIn("def requireModelDirectory", main)
        self.assertIn("mainCheckpointRoot = requireModelDirectory", main)
        self.assertNotIn("mainCheckpointRoot = validateHostDirectory", main)
        self.assertIn("must be a pre-existing directory", main)

    def test_legacy_csv_defaults_to_clinical_notes(self):
        main = (PIPEVAR / "main.nf").read_text()
        self.assertIn(
            "!(params.note != null && clean_note == 'no')",
            main,
        )


if __name__ == "__main__":
    unittest.main()
